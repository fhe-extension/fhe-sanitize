//#pragma once

#include <stdlib.h>
#include "tlwe.h"
//#include "random.h"
//#include "types.h"


/* TGSW samples memory allocation and freeing */
void tgsw_sample_init(tgsw_sample *ct);
void tgsw_sample_clear(tgsw_sample ct);

/* TGSW set to zero */
void tgsw_sample_zero(tgsw_sample *ct);

/* Encrypts polynomial m using key sk and noise parameter param */
tgsw_sample tgsw_encrypt(tlwe_sk sk, double param, int* m);
void tgsw_encrypt_over(tgsw_sample ct, tlwe_sk sk, double param, int* m);

/* Decrypts ciphertext ct using key sk */
int* tgsw_decrypt(tlwe_sk sk, tgsw_sample ct);
void tgsw_decrypt_over_and_keep(int *m, tlwe_sk sk, tgsw_sample ct);

/* Generates a boostrapping key */
bsk generate_bsk(lwe_sk sk_in, tlwe_sk sk_out, double param, int n);
void bsk_clear(bsk bsk, int n);

/* Exeternal product between a TLWE sample and a TGSW sample */
tlwe_sample external_product(tlwe_sample tlwe, tgsw_sample tgsw);

/* Blind_rotate function used during bootstrapping */
tlwe_sample blind_rotate(lwe_sample lwe, int n, uint64_t* test_vector, tgsw_sample *bsk);

/* In place bootstrapping of an LWE ciphertext */
void bootstrap(lwe_sample ct, bsk bsk, ksk ksk, int M_out, int n);
