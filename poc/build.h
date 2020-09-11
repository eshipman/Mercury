#ifndef BUILD_H
#define BUILD_H

#include <stdint.h>

/**
 * FUNCTION: simulate_gsm
 * ----------------------
 * Simulate a GSM compressed voice channel by encoding and decoding the signal
 * using AMR-NB.
 * ARGS  : signal - The signal to run through simulated GSM
 * RETURN:
 * INPUT :
 * OUTPUT: signal - The signal as it would be received on the other end of a
 *                  GSM call over AMR-NB
 */
void simulate_gsm(double *signal);

/**
 * FUNCTION: get_frequency_domain
 * ------------------------------
 * Convert random complex numbers into a properly structured frequency domain
 * that will result in a real-valued time domain.
 * ARGS  : z - The random complex numbers
 *         N - The size of the alphabet
 *         M - The number of samples per symbol
 *         K - The number of active subcarriers
 * RETURN: The frequency domain generated from the complex numbers
 * INPUT :
 * OUTPUT:
 */
complex double* get_frequency_domain(complex double *z, int N, int M, int K);

/**
 * FUNCTION: encode
 * ----------------
 * Encode an integer as a symbol from the alphabet.
 * ARGS  : value - The integer to encode
 *         syms  - The alphabet to encode the integer in
 * RETURN: The corresponding symbol from the alphabet
 * INPUT :
 * OUTPUT:
 */
double* encode(int value, double **syms);

/**
 * FUNCTION: decode
 * ----------------
 * Decode a waveform using the given alphabet.
 * ARGS  : recv - The waveform received
 *         syms - The alphabet to use for decoding
 *         N    - The size of the alphabet
 *         M    - The number of samples per symbol
 * RETURN: The decoded value
 * INPUT :
 * OUTPUT:
 */
int decode(double *recv, double **syms, int N, int M);

/**
 * FUNCTION: state_to_time
 * -----------------------
 * Convert a state of complex numbers to a time domain waveform.
 * ARGS  : state - The state (active fourier bins/subcarriers)
 *         N     - The size of the alphabet
 *         M     - The number of samples per symbol
 *         K     - The number of active subcarriers
 * RETURN: The state represented as a list of real-valued time domain waveforms
 * INPUT :
 * OUTPUT:
 */
double** state_to_time(complex double *state, int N, int M, int K);

/**
 * FUNCTION: optimize
 * ------------------
 * Optimize the initial alphabet to survive compression using Simulated
 * Annealing (SA).
 * ARGS  : x     - The initial alphabet's frequency domain
 *         N     - The number of symbols
 *         K     - The number of active subcarriers
 *         delta - The "closeness" factor (the maximum delta between values of
 *                 a component for neighbor generation)
 *         k_max - The maximum number of iterations
 * RETURN: An alphabet of Speech-Like Symbols optimized to survive compression
 * INPUT :
 * OUTPUT:
 */
double** optimize(complex double **x, int N, int M, int K, int delta, int k_max);

#ifndef RANDOM_NOISE_ITER
#define RANDOM_NOISE_ITER 10000
#endif

#ifndef MIN_N
#define MIN_N 1
#endif

#ifndef MAX_N
#define MAX_N 256
#endif

#ifndef MIN_M
#define MIN_M 1
#endif

#ifndef MAX_M
#define MAX_M 256
#endif

#ifndef MIN_K
#define MIN_K 1
#endif

#ifndef MAX_K
#define MAX_K 128
#endif

#define SUCCESS     0
#define ERROR_USAGE 1

#endif
