#ifndef BUILD_H
#define BUILD_H

#include <stdint.h>

/**
 * STRUCTURE: PSArgs
 * -----------------
 * The required arguments for Pattern Search.
 * D         - The list of directional vectors
 * delta_0   - The initial delta
 * tolerance - The cost at which to stop
 * max_iter  - The maximum number of iterations
 * N         - The number of symbols in the alphabet
 * K         - The number of active subcarriers
 */
struct PSArgs {
    complex double **D;
    double delta_0;
    double tolerance;
    int max_iter;
    int N;
    int M;
    int K;
};

/**
 * FUNCTION: encode
 * ---------------
 * Encode a value as a symbol of the alphabet.
 * ARGS  : value - The integer value to encode
 *         syms  - The alphabet to use for encoding
 * RETURN: The encoded symbol
 * INPUT :
 * OUTPUT:
 */
double*
encode(int value, double **syms);

/**
 * FUNCTION: decode
 * ---------------
 * Decode a symbol to a value using the given alphabet.
 * ARGS  : recv - The received signal
 *         syms - The alphabet to use
 *         N    - The number of symbols in the alphabet
 *         M    - The number of samples per symbol
 * RETURN: The most likely decoded value
 * INPUT :
 * OUTPUT:
 */
int
decode(double *recv, double **syms, int N, int M);

/**
 * FUNCTION: initial_state
 * -----------------------
 * Generate an initial pseudorandom state to be optimized.
 * ARGS  : N - The number of symbols in the alphabet
 *         M - The number of samples per symbol
 *         K - The number of active subcarriers
 * RETURN: A random state vector with dimensionality N*K
 * INPUT :
 * OUTPUT:
 */
complex double*
initial_state(int N, int M, int K);

/**
 * FUNCTION: get_D
 * ---------------
 * Returns a list of directional vectors. The result must be a positve spanning
 * set for R^(N*K). If no assumptions are made about the solution space
 * structure, the pattern D = [I, -I], where I is an identity matrix with sides
 * of 2*N*K, will suffice.
 * ARGS  : N - The number of symbols in the alphabet
 *         K - The number of active subcarriers
 * RETURN: A list of directional vectors (coordinate vectors)
 * INPUT :
 * OUTPUT:
 */
complex double**
get_D(int N, int K);

/**
 * FUNCTION: score
 * ---------------
 * Determine the error rate for the given state by encoding 1000 random values,
 * compression, then decompressing, and finally decoding the symbols.
 * ARGS  : state - The state vector to score
 *         N     - The number of symbols in the alphabet
 *         M     - The number of samples per symbol
 *         K     - The number of active subcarriers
 * RETURN: The error rate in the range [0,1]
 * INPUT :
 * OUTPUT:
 */
double
score(complex double *state, int N, int M, int K);

/**
 * FUNCTION: is_feasible
 * ---------------------
 * Determines if the given directional vector can feasibly be applied to the
 * given state vector using the delta.
 * ARGS  : state - The state vector
 *         size  - The number of elements in the state
 *         d     - The directional vector
 *         delta - The distance to move away from the state for each component
 * RETURN: 0 if infeasible, 1 if feasible
 * INPUT :
 * OUTPUT:
 */
int
is_feasible(complex double *state, int size, complex double *d, double delta);

/**
 * FUNCTION: set_neighbor
 * ----------------------
 * Set the candidate state vector by applying the directional vector and delta
 * to the current state vector.
 * ARGS  : dest  - The destination to store the candidate
 *         state - The current state vector
 *         size  - The number of elements in the state
 *         d     - The directional vector
 *         delta - The distance to move away from the state for each component
 * RETURN:
 * INPUT :
 * OUTPUT:
 */
void
set_neighbor(complex double *dest, complex double *state, int size,
             complex double *d, double delta);

/**
 * FUNCTION: optimize
 * ------------------
 * Optimize the starting state using Pattern Search and the given parameters
 * ARGS  : state - The starting state to optimize
 *         args  - The arguments required for Pattern Search
 * RETURN: A real-valued time domain for each symbol of the alphabet
 * INPUT :
 * OUTPUT:
 */
double**
optimize(complex double *state, struct PSArgs args);

/**
 * FUNCTION: E
 * -----------
 * Determine the error rate of the given state by checking the given number of
 * symbols. A random value is generated from 0..iter. They are all
 * concatenated, compressed and decompressed using AMR-NB, and decoded.
 * ARGS  : state - The state to get the error rate for
 *         N     - The number of symbols in the alphabet
 *         M     - The number of samples per symbol
 *         K     - The number of active subcarriers
 *         iter  - The number of random symbols to encode
 * RETURN: The error rate in the range [0,1]
 * INPUT :
 * OUTPUT:
 */
double E(complex double *state, int N, int M, int K, int iter);

/**
 * FUNCTION: state_to_time
 * -----------------------
 * Convert the state vector to the real-valued time domain representation. The
 * result is the alphabet of waveforms as floating-point PCM.
 * ARGS  : state - The state vector of active fourier bins
 *         N     - The number of symbols in the alphabet
 *         M     - The number of samples per symbol
 *         K     - The number of active subcarriers
 * RETURN: A time domain representation of the alphabet
 * INPUT :
 * OUTPUT:
 */
double** state_to_time(complex double *state, int N, int M, int K);

/**
 * FUNCTION: get_time_domain
 * -------------------------
 * Convert the frequency domain representation of the alphabet to a time domain
 * representation.
 * ARGS  : phi - The frequency domain for each symbol
 *         N   - The number of symbols in the alphabet
 *         M   - The number of samples per symbol
 *         K   - The number of active subcarriers
 * RETURN: The time domain representation of the alphabet
 * INPUT :
 * OUTPUT:
 */
double** get_time_domain(complex double **phi, int N, int M, int K);

/**
 * FUNCTION: to_time_domain
 * ------------------------
 * Convert the frequency domain representation of a single symbol to a time
 * domain representation.
 * ARGS  : phi - The freqency domain representation of the symbol
 *         N   - The number of symbols in the alphabet
 *         M   - The number of samples per symbol
 *         K   - The number of active subcarriers
 * RETURN: The time domain representation of the symobl
 * INPUT :
 * OUTPUT:
 */
double* to_time_domain(complex double *phi, int N, int M, int K);

/**
 * FUNCTION: get_frequency_domain
 * ------------------------------
 * Convert active fourier bins to a frequency domain representation.
 * ARGS  : z - Random complex numbers describing the active fourier bins.
 *         N - The number of symbols in the alphabet
 *         M - The number of samples per symbol
 *         K - The number of active subcarriers
 * RETURN: The frequency domain for a symbol
 * INPUT :
 * OUTPUT:
 */
complex double* get_frequency_domain(complex double *z, int N, int M, int K);

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

#ifndef THREADS
#define THREADS 1
#endif

#endif
