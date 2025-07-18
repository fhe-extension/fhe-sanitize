#pragma once

#include <math.h>
#include "random.h"
#include "types.h"

/* LWE samples memory allocation and freeing */
void lwe_sample_init(lwe_sample *ct, int n);
void lwe_sample_clear(lwe_sample ct);

/* LWE sample set to zero */
void lwe_sample_zero(lwe_sample ct, int n);

/* Binary secret key of size n allocation, generation and freeing */
void lwe_sk_init(lwe_sk *s, int n);
lwe_sk lwe_keygen(int n);
lwe_sk lwe_gaussian_keygen(int n, double param);
void lwe_sk_clear(lwe_sk s);

/* Encrypts message mu under secret key s with message space M and noise parameter param, n is the dimension */
lwe_sample lwe_encrypt(lwe_sk s, int M, double param, int mu, int n);
void lwe_encrypt_over(lwe_sample ct, lwe_sk s, int M, double param, int mu, int n);

/* Decrypts ciphertext ct with secret key s and message space M, n is the dimension */
int lwe_decrypt(lwe_sk s, int M, lwe_sample ct, int n);
int lwe_decrypt_and_keep(lwe_sk s, int M, lwe_sample ct, int n, double *err);

/* Generation and freeing of key switching keys */
ksk generate_ksk(lwe_sk sk_in, lwe_sk sk_out, double param, int n_in, int n_out);
void ksk_clear(ksk ksk, int n_in);

/* Switches from LWE samples of dimension n_in to LWE samples of dimensions n_out */
lwe_sample keyswitch(ksk ksk, int n_in, int n_out, lwe_sample in);
void keyswitch_over_and_keep(lwe_sample out, ksk ksk, int n_in, int n_out, lwe_sample in);
