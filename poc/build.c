#ifndef BUILD_C
#define BUILD_C

#include <string.h>     /* memcpy() */
#include <stdio.h>      /* Printing to stdout/stderr */
#include <stdlib.h>     /* malloc()/free() */
#include <complex.h>    /* Complex numbers */
#include <errno.h>      /* Errors from parsing user input */

#include "build.h"
#include "utils.c"

/* Encode a symbol in the alphabet */
double* encode(int value, double **syms) {
    return syms[value];
}

/* Decode by finding which has the maximum dot product */
int decode(double *recv, double **syms, int N, int M) {
    int i, j;
    double prod, max;
    int output = 0;

    /* Check each symbol */
    for (i = 0; i < N; i++) {
        prod = 0;

        /* Perform the dot product of this symbol with the received symbol */
        for (j = 0; j < M; j++) {
            prod += syms[i][j] * recv[j];
        }

        /* If the dot product is higher, save this index */
        if (i == 0 || prod > max) {
            max = prod;
            output = i;
        }
    }

    /* Return the most likely symbol */
    return output;
}

/* Generate the initial state */
complex double* initial_state(int N, int M, int K)
{
    complex double *state;  /* The generated initial state */
    int i;  /* Loop counters */

    /* Allocate N*K complex numbers for the state */
    state = malloc(N * K * sizeof(complex double));

    /* Generate K active fourier bins for each of the N symbols */
    for (i = 0; i < N * K; i++)
        state[i] = (2 * drandom() - 1) + (2 * drandom() - 1) * I;

    return state;
}

/* Get the directional vectors */
complex double** get_D(int N, int K)
{
    int i, j;   /* Loop counters */
    complex double **output;    /* The directional vectors */

    /* Initialize the output with 2 * P vectors */
    output = (complex double**) malloc(sizeof(complex double*) * 2 * (N * 2 * K));

    for (i = 0; i < 2 * (N * 2 * K); i++) {
        /* Allocate each vector as N*K dimensions */
        output[i] = (complex double*) malloc(N * K * sizeof(complex double));

        /* Set everything to 0 */
        for (j = 0; j < N * K; j++)
            output[i][j] = 0.0;

        /* Set the first half to positive, the second is negative */
        if (i < N * 2 * K) {
            /* Odd positions are real, even are imaginary */
            if (i % 2 == 0)
                output[i][i / 2] = 1.0;
            else
                output[i][i / 2] = 1.0 * I;
        } else {
            /* Odd positions are real, even are imaginary */
            if (i % 2 == 0)
                output[i][(i % (N * 2 * K)) / 2] = -1.0;
            else
                output[i][(i % (N * 2 * K)) / 2] = -1.0 * I;
        }
    }

    return output;
}

/* Score the state */
double score(complex double *state, int N, int M, int K)
{
    return E(state, N, M, K, 1000);
}

/* Determine if the state and directional vector combo is feasible */
int is_feasible(complex double *state, int size, complex double *d, double delta)
{
    complex double tmp; /* Temporary checking variable */
    int fail,   /* Failure checking variable */
        i;      /* Loop counter */

    fail = 0;

    /* Check each component of the state vector */
    for (i = 0; i < size; i++) {
        /* Calculate the new component */
        tmp = state[i] + (d[i] * delta);

        /* If value is outside either bound, this state fails */
        fail = creal(tmp) < -1.0 || creal(tmp) > 1.0
            || cimag(tmp) < -1.0 || cimag(tmp) > 1.0;

        /* If this component failed, return early */
        if (fail)
            return 0;
    }

    /* The state must not have failed */
    return 1;
}

/* Set the next candidate vector by applying the directional vector to the state */
void set_neighbor(complex double *dest, complex double *state, int size, complex double *d, double delta)
{
    int i;  /* Loop counter */

    /* Set every component of the state vector */
    for (i = 0; i < size; i++) {
        /* Set the destination to state + (d * delta) */
        dest[i] = state[i] + (d[i] * delta);
    }
}

/* Optimize the starting state via Pattern Search */
double** optimize(complex double *state, struct PSArgs args, double *ber)
{
    int min_index,  /* The index of the minimum scoring candidate state */
        n,      /* Dimensionality of the state vector */
        P,      /* The number of directional vectors */
        k, i;   /* Loop counter */
    double delta,       /* The delta between iterations */
           cost,        /* The current state's cost score */
           min_cost;    /* The minimum cost score of the candidate states */
    complex double **next;  /* The next iteration's state vectors */

    /* Set the dimensionality */
    n = args.N * args.K * 2;

    /* Set the number of directions */
    P = 2 * n;

    /* Set the initial delta */
    delta = args.delta_0;

    /* Get the initial score */
    cost = score(state, args.N, args.M, args.K);

    /* Allocate space for P state vectors of dimensionality n */
    next = (complex double**) malloc(sizeof(complex double*) * P);
    for (i = 0; i < P; i++) {
        next[i] = (complex double*) malloc(sizeof(complex double) * args.N
                * args.K);
    }

    /* Loop until either the max iterations or target cost is reached */
    for (k = 0; k < args.max_iter && cost >= args.tolerance; k++) {

        min_index = -1;
        min_cost = 1.0;

        /* Check each candidate state vector */
#pragma omp parallel for num_threads(THREADS)
        for (i = 0; i < P; i++) {
            /*
             * If the candidate is feasible, set the next state
             * Else, set it the first component to 2 to mark is infeasible
             */
            if (is_feasible(state, args.N * args.K, args.D[i], delta)) {
                /* Set the state */
                set_neighbor(next[i], state, args.N * args.K, args.D[i], delta);

                /* Get its score */
                double tmp = score(next[i], args.N, args.M, args.K);

                /* If it's the best score so far, save its info */
                if (i == 0 || tmp < min_cost) {
                    min_index = i;
                    min_cost = tmp;
                }
            }
        }

        /*
         * If there was a better state, choose it
         * If not, halve the search distance
         */
        if (min_index < 0) {
            printf("%06d, %.3e: No suitable candidate state vectors\n", k,
                    delta);
            
            delta /= 2.0;
        } else if (min_index >= 0 && min_cost < cost) {
            printf("%06d, %.3e: Changing from %.2e to %.2e\n", k, delta, cost,
                    min_cost);
            
            /* Set the cost */
            cost = min_cost;

            /* Copy each element of the state */
            for (i = 0; i < args.N * args.K; i++) {
                state[i] = next[min_index][i];
            }

            /* Double the search distance */
            delta *= 2.0;
        } else {
            printf("%06d, %.3e: Keeping %.2e instead of %.2e\n", k, delta, cost,
                    min_cost);

            /* Halve the search distance */
            delta /= 2.0;
        }

    }

    /* Store the BER */
    *ber = cost;

    /* Clean up memory */
    for (i = 0; i < P; i++)
        free(next[i]);
    free(next);

    return state_to_time(state, args.N, args.M, args.K);
}

/* Calculate the error rate */
double E(complex double *state, int N, int M, int K, int iter)
{
    double **td,    /* The time domain of the alphabet */
           output,  /* The error rate score */
           *tmp;    /* Temporary */
    int *values,    /* The list of values encoded */
        wrong,      /* The number of incorrect decodes */
        total,      /* The total number of decodes */
        i, j;       /* Loop counters */

    output = 0.0;

    /* Get the time domain representation of the state vector */
    td = state_to_time(state, N, M, K);

    /* Allocate space for the values and their encodings */
    tmp = malloc(iter * M * sizeof(double));
    values = malloc(iter * sizeof(int));

    /* Encode iter random symbols */
    for (i = 0; i < iter; i++) {
        /* Pick a random value to encode */
        values[i] = lrandom() % N;

        /* Copy the value from the alphabet to the signal array */
        memcpy(&tmp[i * M], td[values[i]], M * sizeof(double));
    }

    /* Compress and decompress the signal to simulate a GSM network */
    simulate_gsm(tmp, iter * M);

    wrong = total = 0;

    /* Try decoding every symbol */
    for (i = 0; i < iter; i++) {
        /* Decode it as an integer */
        int dec = decode(&tmp[i * M], td, N, M);

        /* Check each of the bits in the decoded value */
        for (j = 0; j < ceil(log(N) / log(2)); j++) {
            /* If the bits don't match, it's wrong */
            if (((dec >> i) & 1) != ((values[i] >> i) & 1))
                wrong++;

            total++;
        }
    }

    /* Return the result as the percentage incorrect */
    output = (double) wrong / total;

    /* Deallocate the no longer needed memory */
    for (i = 0; i < N; i++)
        free(td[i]);
    free(tmp);
    free(td);

    return output;
}

/* Convert a state vector to a time domain representation of the alphabet */
double** state_to_time(complex double *state, int N, int M, int K)
{
    complex double **fd;    /* The frequency domain of the state */
    int i, j;   /* Loop counters */

    /* Allocate memory for N symbols */
    fd = (complex double**) malloc(N * sizeof(complex double));

    /* Convert each of the K subcarriers into frequency domains */
    for (i = 0; i < N; i++)
        fd[i] = get_frequency_domain(&state[i * K], N, M, K);

    /* Convert the frequency domains to time domains */
    return get_time_domain(fd, N, M, K);
}

/* Convert an alphabet from frequency domain to time domain */
double** get_time_domain(complex double **phi, int N, int M, int K)
{
    double **output;    /* The real-valued time domain */
    int i;   /* Loop counter */

    /* Allocate memory for the real-valued time domain */
    output = (double**) malloc(N * sizeof(double*));

    /* Convert each symbol to its time domain */
    for (i = 0; i < N; i++)
        output[i] = to_time_domain(phi[i], N, M, K);

    return output;
}

/* Convert a frequency domain to a real-valued time domain */
double* to_time_domain(complex double *phi, int N, int M, int K)
{
    complex double *g;  /* The time domain */
    double *output;     /* The real-valued time domain */
    int k;   /* Loop counters */

    /* Allocate memory for the output of M samples */
    output = (double*) malloc(M * sizeof(double));

    /* Get the time domain via the inverse DFT */
    g = inverse_dft(phi, M + 1);

    /* Normalize the symbol power */
    for (k = 0; k < M + 1; k++) {
        /* If one component is non-zero, divide by norm() */
        if (creal(g[k]) != 0 || cimag(g[k]) != 0)
            g[k] = g[k] / norm(g[k]);

        /* Save the real part */
        if (k < M)
            output[k] = creal(g[k]);
    }

    /* De-allocate the complex time domain */
    free(g);

    return output;
}

/* Convert random complex numbers to a frequency domain */
complex double* get_frequency_domain(complex double *z, int N, int M, int K)
{
    complex double *phi;    /* The frequency domain */
    int k_N,    /* The index corresponding to the Nyquist frequency */
        k;      /* Loop counter */

    /* The index of the Nyquist frequency is M / 2 */
    k_N = M / 2;

    /* Allocate memory for the frequency domain */
    phi = (complex double*) malloc((M + 1) * sizeof(complex double));

    /*
     * Produce the complex spectrum for the symbol according to the following:
     *
     * phi_k =  0,          if k = 0                    (1)
     *          z_(k-1),    if k = 1, ..., K            (2)
     *          0,          if k = K + 1, ..., k_N      (3)
     *          z_(M-k-1)*  if k = k_N + 1, ..., M - 1  (4)
     *          0           if k = M
     */

    /* Rule 1: phi_k = 0, if k = 0 */
    phi[0] = 0;

    /* Rule 2: phi_k = z_(k-1), if k = 1, ..., K */
    for (k = 1; k <= K; k++)
        phi[k] = z[k - 1];

    /* Rule 3: phi_k = 0, if k = K + 1, ..., k_N */
    for (k = K + 1; k <= k_N; k++)
        phi[k] = 0;

    /* Rule 4: phi_k = z_(M-k-1)* if k = k_N + 1, ..., M - 1 */
    for (k = k_N + 1; k < M; k++)
        phi[k] = conj(z[M - k - 1]);

    phi[M] = 0;

    return phi;
}

int main(int argc, char **argv)
{
    int N,      /* The number of symbols in the alphabet */
        M,      /* The number of samples per symbol */
        K;      /* The number of subcarriers */
    char *ptr;  /* Pointer for passing to strtol() */
    double **alphabet;  /* The speech like symbols */
    double ber; /* The BER for this alphabet */

    /* Verify the number of arguments is correct */
    if (argc < 4) {
        fprintf(stderr,
                "Usage: %s <N value ([1,inf])> <M value ([1,inf])> <K value ([1,inf])>\n",
                argv[0]);
        return ERROR_USAGE;
    }

    /* Parse the N-value, handling errors */
    ptr = 0;
    N = strtol(argv[1], NULL, 10);

    /* If errors were encountered, print the error and exit */
    if (ptr == argv[1] || errno == ERANGE || errno == EINVAL) {
        perror("strtol");
        return ERROR_USAGE;
    } else if (N < MIN_N || N > MAX_N) {
        fprintf(stderr, "Error: N value (%d) is out of bounds [%d,%d]\n", N,
                MIN_N, MAX_N);
        return ERROR_USAGE;
    }

    /* Parse the M-value, handling errors */
    ptr = 0;
    M = strtol(argv[2], NULL, 10);

    /* If errors were encountered, print the error and exit */
    if (ptr == argv[2] || errno == ERANGE || errno == EINVAL) {
        perror("strtol");
        return ERROR_USAGE;
    } else if (M < MIN_M || M > MAX_M) {
        fprintf(stderr, "Error: M value (%d) is out of bounds [%d,%d]\n", M,
                MIN_M, MAX_M);
        return ERROR_USAGE;
    }

    /* Parse the K-value, handling errors */
    ptr = 0;
    K = strtol(argv[3], NULL, 10);

    /* If errors were encountered, print the error and exit */
    if (ptr == argv[3] || errno == ERANGE || errno == EINVAL) {
        perror("strtol");
        return ERROR_USAGE;
    } else if (K < MIN_K || K > MAX_K) {
        fprintf(stderr, "Error: K value (%d) is out of bounds [%d,%d]\n", K,
                MIN_K, MAX_K);
        return ERROR_USAGE;
    } else if (K >= M / 2) {
        fprintf(stderr, "Error: K value must be less than M/2\n");
        return ERROR_USAGE;
    }

    /* Build the arguments */
    struct PSArgs args;
    args.D = get_D(N, K);   /* Get the standard directional vectors */
    args.delta_0 = 0.5;    /* Set the initial delta */
    args.tolerance = 0.0001; /* Set the tolerance */
    args.max_iter = 10000;   /* Set the maximum number of iterations */
    args.N = N;     /* Set the given N, M, & K */
    args.M = M;
    args.K = K;

    /* Generate an initial state */
    complex double *start = initial_state(N, M, K);

    /* Optimize the state and return a real-valued time domain alphabet */
    alphabet = optimize(start, args, &ber);

    printf("-----BEGIN-HEADER-----\n");
    printf("#ifndef SYMBOLS_H\n");
    printf("#define SYMBOLS_H\n");
    printf("\n\n/*\n");
    printf("  These symbols have a BER of %e%%\n", ber * 100.0);
    printf("\n*/\n\n");
    printf("#define N %d\n", N);
    printf("#define N %d\n\n", M);
    printf("double SYMBOLS[N][M] = {\n");
    for (int i = 0; i < N; i++) {
        printf("    {%lf", alphabet[i][0]);
        for (int j = 1; j < M; j++) {
            printf(", %lf", alphabet[i][j]);
        }
        printf("}");
        if (i < N - 1)
            printf(",\n");
    }
    printf("};\n\n");
    printf("#endif\n");
    printf("-----END-HEADER-----");

    return SUCCESS;
}
#endif
