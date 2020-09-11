#ifndef BUILD_C
#define BUILD_C

#include <string.h>     /* memcpy() */
#include <stdio.h>      /* Printing to stdout/stderr */
#include <stdlib.h>     /* malloc()/free() */
#include <complex.h>    /* Complex numbers */
#include <errno.h>      /* Errors from parsing user input */

#include <opencore-amrnb/interf_enc.h>
#include <opencore-amrnb/interf_dec.h>

#include "build.h"
#include "utils.c"      /* randint() */

void simulate_gsm(double *signal)
{
    /*
    enum Mode mode = MR122;
    void *amr;
    uint8_t* outbuff;
    short *inbuff;
    int16_t *tmp;
    double *output;
    uint8_t* decbuff;
    int dtx = 0;
    int i;

    inbuff = malloc(length * sizeof(short));
    for (i = 0; i < length; i++) {
        inbuff[i] = signal[i] * 32767;
    }
    
    outbuff = malloc(length);

    amr = Encoder_Interface_init(dtx);

    for (i = 0; i < length / 160; i++) {
        Encoder_Interface_Encode(amr, mode, inbuff, outbuff, 0);
    }

    Encoder_Interface_exit(amr);

    amr = Decoder_Interface_init();

    Decoder_Interface_Decode(amr, outbuff, tmp, 0);

    free(inbuff);
*/
    return;
}

complex double** generate_alphabet(int N, int M, int K)
{
    double complex *z,      /* The complex numbers for generating the symbols */
                   **phi;   /* The frequency domain for each symbol */
    int i, k;   /* Loop counters */

    /* Allocate the complex numbers */
    z = (double complex*) malloc(K * sizeof(double complex));

    /* Allocate the frequency domains (an array for each symbol) */
    phi = (double complex**) malloc(N * sizeof(double complex*));

    /* Generate N symbols */
    for (i = 0; i < N; i++) {

        /* Generate K random complex numbers for each symbol */
        for (k = 0; k < K; k++) {
            z[k] = (2 * drandom() - 1) + (2 * drandom() - 1) * I;
        }

        phi[i] = get_frequency_domain(z, N, M, K);
    }

    free(z);

    return phi;
}

/* Generate a neighbor to the given state vector */
complex double* neighbor(complex double *current, int n, double delta)
{
    complex double *output, /* The neighboring state vector */
                   tmp;     /* A temp variable to store the possible element */
    int done,   /* Loop condition */
        i;      /* Loop counter */

    /* Allocate the neighbor's state vector */
    output = (complex double*) malloc(n * sizeof(complex double));

    /* Set each element in the to be near each element of the curretn */
    for (i = 0; i < n; i++) {
        done = 0;

        /* Generate random values until one is found within the bounds */
        while (!done) {
            /*
             * Add the random complex number to the current state. The
             * components must be within delta of the current value.
             */
            tmp = current[i]
                + (complex double) ((delta * drandom() * 2 - 1)
                + (delta * drandom() * 2 - 1) * I);

            /*
             * This search is done when both the real and imaginary components
             * are in the range [-1, 1].
             */
            done = creal(tmp) >= -1 && creal(tmp) <= 1
                && cimag(tmp) >= -1 && cimag(tmp) <= 1;
        }

        printf("neighbor: <%lf, %lf>\n", creal(tmp - current[i]), cimag(tmp - current[i]));
        /* Set the output value for this index */
        output[i] = tmp;
    }

    return output;
}

/*
 * Get the probability of choosing a solution based on the error rates and time
 */
double P(double old, double new, double t)
{
    return exp(-new / (old * t));
}

/* Convert random complex numbers to a frequency domain */
complex double* get_frequency_domain(complex double *z, int N, int M, int K)
{
    complex double *phi;    /* The frequency domain */
    int k;      /* Loop counter */
    int k_N;    /* The index corresponding to the Nyquist frequency */

    /* Is k_N always M / 2 ? */
    k_N = M / 2;

    /* Allocate space for the frequency domain */
    phi = (complex double*) malloc((M + 1) * sizeof(complex double));

    /*
     * Produce the complex spectrum for the symbol according to the following:
     *
     * phi_k =  0,          if k = 0                (1)
     *          z_k,        if k = 1, ..., K        (2)
     *          0,          if k = K + 1, ..., k_N  (3)
     *          z_(M-k+1])* if k = k_N + 1, ..., M  (4)
     */

    /* Rule 1: phi_k = 0, if k = 0 */
    phi[0] = 0;

    /* Rule 2: phi_k = z_k, if k = 1, ..., K */
    for (k = 1; k <= K; k++) {
        phi[k] = z[k];
    }

    /* Rule 3: phi_k = 0, if k = K + 1, ..., k_N */
    for (k = K + 1; k <= k_N; k++) {
        phi[k] = 0;
    }

    /* Rule 4: phi_k = z_(M-k+1)* if k = k_N + 1, ..., M */
    for (k = k_N + 1; k <= M; k++) {
        phi[k] = conj(z[M - k + 1]);
    }

    return phi;
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

double** get_time_domain(complex double **phi, int N, int M, int K)
{
    double **output;    /* The real-valued time domain */
    int i;   /* Loop counter */

    /* Allocate space for the real-valued time domain */
    output = (double**) malloc(N * sizeof(double*));

    /* Convert each symbol to its time domain */
    for (i = 0; i < N; i++) {
        output[i] = to_time_domain(phi[i], N, M, K);
    }

    return output;
}

double E(complex double *state, int N, int M, int K, int iter)
{
    double **td;
    double output,
           *tmp;
    int *values;
    int i;
    int wrong, total;

    output = 0.0;

    td = state_to_time(state, N, M, K);

    tmp = malloc(iter * M * sizeof(double));
    values = malloc(iter * sizeof(int));

    for (i = 0; i < iter; i++) {
        /* Pick a random value to encode */
        values[i] = lrandom() % N;

        /* Copy the value from the alphabet to the signal array */
        memcpy(&tmp[i * M], td[values[i]], M * sizeof(double));
    }

    /* Simulate a GSM network */
    simulate_gsm(tmp);

    wrong = total = 0;

    for (i = 0; i < M; i++) {
        printf("E: %lf\n", tmp[i]);
    }

    for (i = 0; i < iter; i++) {
        int dec = decode(&tmp[i * M], td, N, M);
        if (dec != values[i]) {
            printf("E: Thought %d, was %d\n", dec, values[i]);
            wrong++;
        }
        total++;
    }

    output = (double) wrong / total;

    for (i = 0; i < N; i++) {
        free(td[i]);
    }
    free(tmp);
    free(td);

    return output;
}

/* Convert a frequency domain to state */
void set_state(complex double *state, complex double **x, int N, int K)
{
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < K; j++) {
            state[i * K + j] = x[i][j];
            printf("set_state: <%lf, %lf>\n", creal(x[i][j]), cimag(x[i][j]));
        }
    }
}

double** state_to_time(complex double *state, int N, int M, int K)
{
    complex double **fd;
    int i;

    fd = (complex double**) malloc(N * sizeof(complex double));

    for (i = 0; i < N; i++)
        fd[i] = get_frequency_domain(state, N, M, K);

    return get_time_domain(fd, N, M, K);
}

double** optimize(complex double **x, int N, int M, int K, int delta, int k_max)
{
    complex double *state,  /* The current state vector */
            *next;   /* The candidate next state */
    double old_ber, /* The score of the current state */
           new_ber; /* The score of the candidate state */
    int ber_valid,  /* Whether the current BER is valid */
        k;          /* Loop counter */

    /* Allocate space for the vector */
    state = (complex double*) malloc(N * K * sizeof(complex double));

    // Set the state
    set_state(state, x, N, K);

    /* Initialize the set flag to 0 */
    ber_valid = 0;

    for (k = 0; k < k_max; k++) {
        /* Generate a candidate neighbor state */
        next = neighbor(state, N * K, delta);

        /* Get the BER of the new state */
        new_ber = E(next, N, M, K, RANDOM_NOISE_ITER);

        /* If the BER isn't valid, update the BER */
        if (!ber_valid)
            old_ber = E(state, N, M, K, RANDOM_NOISE_ITER);

        /* Set to 1 so BER isn't calculated redundantly */
        ber_valid = 1;

        /*
         * Probabilistically choose the generated state based on the time and
         * the BER for both states.
         */
        if (P(old_ber, new_ber, (k + 1) / k_max) >= drandom()) {
            printf("Changing state from %lf to %lf BER\n", old_ber, new_ber);
            memcpy(state, next, N * K * sizeof(complex double));
            ber_valid = 0;
        } else {
            printf("Keeping state of %lf instead of %lf\n", old_ber, new_ber);
        }
        fflush(stdout);

        free(next);
    }

    return state_to_time(state, N, M, K);
}

double* encode(int value, double **syms) {
    return syms[value];
}

/* Decode by finding which has the maximum dot product */
int decode(double *recv, double **syms, int N, int M) {
    int i, j;
    double prod, max;
    int output = 0;
    for (i = 0; i < N; i++) {
        prod = 0;
        for (j = 0; j < M; j++) {
            prod += syms[i][j] * recv[j];
        }
        if (i == 0 || (i > 0 && prod > max)) {
            max = prod;
            output = i;
        }
    }

    return output;
}

int main(int argc, char **argv)
{
    int N,      /* The number of symbols in the alphabet */
        M,      /* The number of samples per symbol */
        K;      /* The number of subcarriers */
    char *ptr;  /* Pointer for passing to strtol() */
    complex double **syms;  /* The speech like symbols */
    double **alphabet;

    /* Verify the number of arguments is correct */
    if (argc < 4) {
        fprintf(stderr,
                "Usage: %s <N value ([1,inf])> <M value ([1,inf])> <K value ([1,inf])>", argv[0]);
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

    syms = generate_alphabet(N, M, K);
    //alphabet = optimize(syms, N, M, K, 0.1, 1000);
    alphabet = get_time_domain(syms, N, M, K);
    
    int i;
    int wrong, total;
    wrong = total = 0;
    for (i = 0; i < 1; i++) {
        int value = lrandom() % N;

        double *sig = encode(value, alphabet);
        int dec = decode(sig, alphabet, N, M);

        if (value != decode(sig, alphabet, N, M)) {
            wrong++;
            printf("Thought %d, was %d\n", dec, value);
        }
        total++;
    }

    printf("Final total %lf%%\n", (double) wrong / total);

    return SUCCESS;
}
#endif
