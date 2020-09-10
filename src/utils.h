#ifndef UTILS_H
#define UTILS_H

#include <stdint.h>

/* Whether or not the PRNG has been seeded yet */
static int PRNG_SEEDED = 0;

/* The internal state of the xorshift* PRNG */
static struct PRNG_STATE {
    uint64_t a;
} PRNG_STATE;

/**
 * FUNCTION: init_prng
 * -------------------
 * Seeds the PRNG with the current time in nanoseconds
 * ARGS  :
 * RETURN:
 * INPUT : The current time in nanoseconds
 * OUTPUT: Updates PRNG_STATE with a 64-bit seed
 */
static void init_prng();

/**
 * Function: lrandom
 * -----------------
 * Generates a pseudorandom number using the xorshift* PRNG
 * ARGS  :
 * RETURN: A pseudorandom 64-bit number
 * INPUT : PRNG_STATE - the current state of the PRNG
 * OUTPUT:
 */
uint64_t lrandom();

/**
 * FUNCTION: drandom
 * -----------------
 * Generates a pseudorandom double in the range [0,1] using the xorshift* PRNG
 * ARGS  :
 * RETURN: A pseudorandom double
 * INPUT : PRNG_STATE - the current state of the PRNG
 * OUTPUT:
 */
double drandom();

/**
 * FUNCTION: dft
 * -------------
 * Compute the Discrete Fourier Transform (DFT) of a set of complex numbers.
 * Algorithm source: https://en.wikipedia.org/wiki/Discrete_Fourier_transform
 * ARGS  : X - The destination to store the result
 *         x - The set of numbers to transform
 *         N - The size of the set
 * RETURN: The transformed complex numbers
 * INPUT :
 * OUTPUT:
 */
complex double* dft(complex double *X, complex double *x, int N);

/**
 * FUNCTION: inverse_dft
 * ---------------------
 * Compute the inverse Discrete Fourier Transform (DFT) of a set of complex
 * numbers.
 * Algorithm source: https://en.wikipedia.org/wiki/Discrete_Fourier_transform
 * ARGS  : x - The destination to store the result
 *         X - The set of numbers to transform
 *         N - The size of the set
 * RETURN: The transformed complex numbers
 * INPUT :
 * OUTPUT:
 */
complex double* inverse_dft(complex double *x, complex double *X, int N);

/**
 * FUNCTION: norm
 * --------------
 * Compute the Euclidean norm of a complex number. The result will always be a
 * real number.
 * ARGS  : z - The complex number
 * RETURN: The Euclidean norm
 * INPUT :
 * OUTPUT:
 */
double norm(complex double z);

#endif
