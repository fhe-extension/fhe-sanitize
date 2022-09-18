#pragma once

//#include <stdlib.h>
//#include <stdio.h>	
//#include <string.h>

//#include "polynomial.h"
//#include "types.h"
#include "lwe.h"
//#include "random.h"
#define _mask ((int64_t) -((uint64_t) 1<<(64-_logBg*_ell)))


/* TLWE samples memory allocation and freeing */
void tlwe_sample_init(tlwe_sample *ct);
void tlwe_sample_clear(tlwe_sample ct);

/* TLWE set to zero */
void tlwe_sample_zero(tlwe_sample ct);

/* TLWE binary secret keys memory allocation, generation and freeing */
void tlwe_sk_init(tlwe_sk *sk);
tlwe_sk tlwe_keygen();
void tlwe_sk_clear(tlwe_sk sk);

/* Adds the tlwe sample a to the tlwe sample out */
void tlwe_add_to(tlwe_sample out, tlwe_sample a);
/* Polynomials multiplications modulo X^N+1 */
void multiply_accumulate_poly(uint64_t *out, uint64_t *a, uint64_t *b);
/* Deterministic decomposition */
void decompose_tlwe(uint64_t *out, tlwe_sample in);

/* Encrypts polynomial m using key sk, message space M and noise parameter/stdev param */
tlwe_sample tlwe_encrypt(tlwe_sk sk, long double param, int M, int* m);
void tlwe_encrypt_over(tlwe_sample ct, tlwe_sk sk, long double param, int M, int* m);
tlwe_sample tlwe_encrypt_zero(tlwe_sk sk, long double param);
void tlwe_encrypt_zero_over(tlwe_sample ct, tlwe_sk sk, long double param);
tlwe_sample tlwe_encrypt_zero_gaussian(tlwe_sk sk, long double param);
void tlwe_encrypt_zero_gaussian_over(tlwe_sample ct, tlwe_sk sk, long double param);

/* Decrypts ciphertext ct using key sk, message space M */
int* tlwe_decrypt(tlwe_sk sk, int M, tlwe_sample ct);
void tlwe_decrypt_over_and_keep(int *m, tlwe_sk sk, int M, tlwe_sample ct);

/* Extracts an LWE secret key from a TLWE secret key */
lwe_sk tlwe_key_extract(tlwe_sk in);
/* Extracts an LWE ciphertext from a TLWE ciphertext */
lwe_sample tlwe_extract(tlwe_sample in);
void tlwe_extract_over_and_keep(lwe_sample out, tlwe_sample in);

/* Generates public keys necessary for sanitization of ciphertexts */
pkg sanitize_pk_gen(tlwe_sk tsk, long double param);
void sanitize_pk_clear(pkg PK);

/* Generates a well distributed TLWE encryption of 0 and adds it to the ciphertext */
void sanitize_pk_enc(tlwe_sample ct, pkg PK);

/* Initialize and generates a pool of encryptions of 0 for use in sanitization */
void precompute_pkenc(pkg PK);
void clear_pkenc();
