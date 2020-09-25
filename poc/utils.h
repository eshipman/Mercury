#ifndef UTILS_H
#define UTILS_H

#include <stdint.h>

/* Whether or not the PRNG has been seeded yet */
static int PRNG_SEEDED = 0;

/* The internal state of the xorshift* PRNG */
/**
 * STRUCTURE: PRNG_STATE
 * ---------------------
 * The current state of the PseudoRandom Number Generator.
 * a - The state as a 64-bit integer
 */
static struct PRNG_STATE {
    uint64_t a;
} PRNG_STATE;

/**
 * FUNCTION: init_prng
 * -------------------
 * Seeds the PRNG with the current time in nanoseconds.
 * ARGS  :
 * RETURN:
 * INPUT : The current time in nanoseconds
 * OUTPUT: Updates PRNG_STATE with a 64-bit seed
 */
static void init_prng();

/**
 * Function: lrandom
 * -----------------
 * Generates a random number using the xorshift* PRNG.
 * ARGS  :
 * RETURN: A pseudorandom 64-bit number
 * INPUT : PRNG_STATE - the current state of the PRNG
 * OUTPUT:
 */
uint64_t lrandom();

/**
 * FUNCTION: drandom
 * -----------------
 * Generates a random double in the range [0,1] using the xorshift* PRNG.
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
 * ARGS  : x - The set of numbers to transform
 *         N - The size of the set
 * RETURN: The transformed complex numbers
 * INPUT :
 * OUTPUT:
 */
complex double* dft(complex double *x, int N);

/**
 * FUNCTION: inverse_dft
 * ---------------------
 * Compute the inverse Discrete Fourier Transform (DFT) of a set of complex
 * numbers.
 * Algorithm source: https://en.wikipedia.org/wiki/Discrete_Fourier_transform
 * ARGS  : X - The set of numbers to transform
 *         N - The size of the set
 * RETURN: The transformed complex numbers
 * INPUT :
 * OUTPUT:
 */
complex double* inverse_dft(complex double *X, int N);

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

/**
 * FUNCTION: simulate_gsm
 * ----------------------
 * Simulate a GSM compressed voice channel by encoding and decoding the signal
 * using AMR-NB at 12.2 Kbps.
 * ARGS  : signal - The signal to run through simulated GSM
 *         length - The number of samples in the signal
 * RETURN:
 * INPUT :
 * OUTPUT: signal - The signal as it would be received on the other end of a
 *                  GSM call over AMR-NB
 */
void simulate_gsm(double *signal, int length);

/**
 * FUNCTION: secded_encode
 * -----------------------
 * Encode the given bytes using SECDED (Single Error Correction, Double Error
 * Detection). The result will always be twice the size of the input.
 * ARGS  : out  - The destination to store the encoded result
 *         data - The data to be encoded
 *         sz   - The size of the input
 * RETURN:
 * INPUT :
 * OUTPUT:
 */
void secded_encode(uint8_t *out, uint8_t *data, size_t sz);

/**
 * FUNCTION: secded_decode
 * -----------------------
 * Decode the given bytes as SECDED (Single Error Correction, Double Error
 * Detection). The result will always be halve the size of the input.
 * ARGS  : out   - The destination to store the decoded data
 *         data  - The data to be decoded
 *         sz    - The size of the input
 *         *cerr - The destination to store the number of corrected errors
 *         *derr - The destination to store the number of detected errors
 * RETURN:
 * INPUT :
 * OUTPUT: *cerr - The number of errors in the input that were corrected
 *         *derr - The number of errors in the input that could not be
 *                 corrected but could be detected
 */
void secded_decode(uint8_t* out, uint8_t *data, size_t sz, int *cerr, int *derr);

#endif
