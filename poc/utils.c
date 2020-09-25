#ifndef UTILS_C
#define UTILS_C

#include <time.h>       /* clock_gettime() for seeding the PRNG */
#include <complex.h>    /* Complex numbers for DFT and its inverse */
#include <math.h>       /* exp() */
#include <stdlib.h>     /* malloc() */

#include <opencore-amrnb/interf_enc.h>  /* For AMR-NB Encoding */
#include <opencore-amrnb/interf_dec.h>  /* For AMR-NB Decoding */

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

/* Simulate a GSM network */
void simulate_gsm(double *signal, int length)
{
    void *amr_enc,  /* The encoding state */
         *amr_dec;  /* The decoding state */
    enum Mode mode = MR122; /* Only encode in 12.2kbps mode for now */
    short enc_in[160];  /* Encoder requires 160 shorts at a time */
    uint8_t enc_out[500];   /* Encoder produces 500 bytes at a time */
    int16_t dec_out[160];   /* Decoder produces 160 16-bit ints ata time */
    int dtx = 1;    /* Simulate with DTX on */

    /* Initialize both states */
    amr_enc = Encoder_Interface_init(dtx);
    amr_dec = Decoder_Interface_init();

    /* Process the input 160 samples at a time */
    for (int i = 0; i < ceil(length / 160); i++) {

        /*
         * Convert to signed 16-bit PCM, padding with zeroes if the signal is
         * less than a multiple of 160 samples.
         */
        for (int j = 0; j < 160; j++) {
            if (i * 160 + j < length)
                enc_in[j] = round(signal[i * 160 + j] * 32767);
            else
                enc_in[j] = 0;
        }

        /* Encode the data, then decode it */
        Encoder_Interface_Encode(amr_enc, mode, enc_in, enc_out, 0);
        Decoder_Interface_Decode(amr_dec, enc_out, dec_out, 0);

        /* Replace the input data with the simulated signal */
        for (int j = 0; j < 160; j++) {
            if (i * 160 + j < length)
                signal[i * 160 + j] = (double) dec_out[j] / 32768.0;
        }
    }
}

/* Encode n bytes to 2n bytes with SECDED ECC */
void secded_encode(uint8_t *out, uint8_t *data, size_t sz)
{
    uint8_t bits[8];    /* The list of the bits for each byte */
    int i, j;   /* Loop counters */
    
    for (i = 0; i < sz; i += 2) {
        uint8_t l8 = 0,
                r8 = 0;
        
        /* Extract each bit */
        for (j = 0; j < 8; j++)
            bits[j] = (data[i] >> (7 - j)) & 0x01;
        
        /* Encode the most significant nibble as the first byte */
        out[i] = (bits[0] ^ bits[1] ^ bits[3]) << 7
               | (bits[0] ^ bits[2] ^ bits[3]) << 6
               | (bits[0])                     << 5
               | (bits[1] ^ bits[2] ^ bits[3]) << 4
               | (bits[1])                     << 3
               | (bits[2])                     << 2
               | (bits[3])                     << 1;

        /* Calculate the additional parity bit */
        l8 =    (out[i] >> 7) ^ (out[i] >> 6) ^ (out[i] >> 5) ^ (out[i] >> 4)
              ^ (out[i] >> 3) ^ (out[i] >> 2) ^ (out[i] >> 1) ^ out[i];

        /* Set the last bit as the additional parity bit */
        out[i] |= (l8 & 0x01);

        /* Encode the least significant nibble as the second byte */
        out[i+1] = (bits[4] ^ bits[5] ^ bits[7]) << 7
                 | (bits[4] ^ bits[6] ^ bits[7]) << 6
                 | (bits[4])                     << 5
                 | (bits[5] ^ bits[6] ^ bits[7]) << 4
                 | (bits[5])                     << 3
                 | (bits[6])                     << 2
                 | (bits[7])                     << 1;;

        /* Calculate the additional parity bit */
        r8 = (out[i+1] >> 7) ^ (out[i+1] >> 6) ^ (out[i+1] >> 5) ^ (out[i+1] >> 4)
           ^ (out[i+1] >> 3) ^ (out[i+1] >> 2) ^ (out[i+1] >> 1) ^ out[i+1];

        /* Set the last bit as the additional parity bit */
        out[i+1] |= (r8 & 0x01);
    }
}

/* Decode 2n bytes of SECDED ECC encoded data into n bytes of data */
void secded_decode(uint8_t* out, uint8_t *data, size_t sz, int *cerr, int *derr)
{
    uint8_t bits[16];   /* The list of bits for each pair */
    int i, j;   /* Loop counters */

    /* Initialize the number of errors found to 0 */
    *cerr = 0;
    *derr = 0;

    /* Loop through each pair of bytes */
    for (i = 0; i < sz; i += 2) {
        uint8_t l8 = 0,
                r8 = 0;
        uint8_t c[2][3],
                c0 = 0,
                c1 = 0;
            
        /* Get each bit individually and get the additional parity */
        for (j = 0; j < 8; j++) {
            bits[j] = (data[i] >> (7 - j)) & 0x01;
            bits[j + 8] = (data[i + 1] >> (7 - j)) & 0x01;
            l8 ^= bits[j];
            r8 ^= bits[j + 8];
        }

        /* Determine the positions of any errors */
        c[0][0] = bits[2] ^ bits[4] ^ bits[6] ^ bits[0];
        c[0][1] = bits[2] ^ bits[5] ^ bits[6] ^ bits[1];
        c[0][2] = bits[4] ^ bits[5] ^ bits[6] ^ bits[3];
        c[1][0] = bits[10] ^ bits[12] ^ bits[14] ^ bits[8];
        c[1][1] = bits[10] ^ bits[13] ^ bits[14] ^ bits[9];
        c[1][2] = bits[12] ^ bits[13] ^ bits[14] ^ bits[11];
        c0 = c[0][0] + 2 * c[0][1] + 4 * c[0][2];
        c1 = c[1][0] + 2 * c[1][1] + 4 * c[1][2];

        /* Decode the data bits */
        out[i / 2] = bits[2]  << 7
                   | bits[4]  << 6
                   | bits[5]  << 5
                   | bits[6]  << 4
                   | bits[10] << 3
                   | bits[12] << 2
                   | bits[13] << 1
                   | bits[14];

        /* If errors occurred in the first byte, make the appropriate flips */
        if (c0 != 0) {
            out[i / 2] ^= 0x01 << c0;
            *cerr += 1;
        }
        
        /* If errors occurred in the second byte, make the appropriate flips */
        if (c1 != 0) {
            out[i / 2] ^= 0x01 << c1;
            *cerr += 1;
        }

        /* Determine the expected additional parity bits */
        l8 ^= bits[2] ^ bits[4] ^ bits[5] ^ bits[6] ^ c[0][0] ^ c[0][1] ^ c[0][2];
        r8 ^= bits[10] ^ bits[12] ^ bits[13] ^ bits[14] ^ c[1][0] ^ c[1][1] ^ c[1][2];

        /* If l8 is non-zero, a double error probably occurred */
        if (l8)
            *derr += 1;
        
        /* If r8 is non-zero, a double error probably occurred */
        if (r8)
            *derr += 1;
    }
}

#endif
