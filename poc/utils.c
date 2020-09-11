#ifndef UTILS_C
#define UTILS_C

#include <time.h>       /* clock_gettime() for seeding the PRNG */
#include <complex.h>    /* Complex numbers for DFT and its inverse */
#include <math.h>       /* exp() */
#include <stdlib.h>     /* malloc() */

#include "utils.h"

/* Initialize the PRNG state */
static void init_prng()
{
    uint64_t time;          /* 64-bit seed, current time in nanoseconds */
    struct timespec curr;   /* The current time's timespec */

    /* Get the current time */
    clock_gettime(CLOCK_MONOTONIC, &curr);

    /* Convert the timespec to the time in nanoseconds */
    time = curr.tv_sec * 10e9 + curr.tv_nsec;

    /* Set the PRNG state to the nanosecond time */
    PRNG_STATE.a = time;
}

/* Generates a random 64-bit number */
uint64_t lrandom()
{
    uint64_t x;

    /* If the PRNG hasn't been seeded, do it */
    if (PRNG_SEEDED == 0) {
        init_prng();
        PRNG_SEEDED = 1;
    }

    /* Apply the xorshift operations */
    x = PRNG_STATE.a;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    PRNG_STATE.a = x;

    /* Return the state multiplied with a constant */
    return x * UINT64_C(0x2545F4914F6CDD1D);
}

/* Generates a random double between 0 and 1 */
double drandom()
{
    /* Generate a random 64-bit number and convert to a double in [0,1] */
    return (double) lrandom() / UINT64_MAX;
}

/* Compute the DFT of a set of complex numbers */
complex double* dft(complex double *x, int N)
{
    complex double *X;  /* The resulting transformed complex numbers */
    int k, n;  /* Loop counters */

    /* Allocate memory for the transformed values */
    X = (complex double*) malloc(N * sizeof(complex double));

    /* Compute the DFT */
    /* Algorithm used from wikipedia */
    for (k = 0; k < N; k++) {
        X[k] = 0;
        for (n = 0; n < N; n++)
            X[k] += x[n] * cexp((-I * 2 * M_PI) / N * k * n);
    }

    /* Return the result */
    return X;
}

/* Compute the inverse DFT of a set of complex numbers */
complex double* inverse_dft(complex double *X, int N)
{
    complex double *x;  /* The resulting transformed complex numbers */
    int k, n;   /* Loop counters */

    /* Allocate memory for the transformed numbers */
    x = (complex double*) malloc(N * sizeof(complex double));

    /* Compute the inverse DFT */
    /* Algorithm used from wikipedia */
    for (n = 0; n < N; n++) {
        x[n] = 0;
        for (k = 0; k < N; k++) {
            x[n] += X[k] * cexp(I * 2 * M_PI * k * n / N);
        }
        x[n] /= N;
    }

    /* Return the result */
    return x;
}

/* Calculate the Euclidean norm of a complex number */
double norm(complex double z)
{
    /* Calculate and return the norm in one line */
    return sqrt(pow(creal(z), 2) + pow(cimag(z), 2));
}

#endif
