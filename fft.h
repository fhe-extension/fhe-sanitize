#ifndef FFT_H
#define FFT_H
#include <complex.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <fftw3.h>
#include "tgsw.h"

typedef fftw_complex* tlwe_sk_fft;
typedef fftw_complex* tgsw_sample_fft;
typedef fftw_complex* pkc_fft;
typedef tgsw_sample_fft *bsk_fft;

/* Debug access only */
extern double *fft_real;
extern fftw_complex *fft_complex;
extern fftw_complex *fft_complex_mid;
extern fftw_complex *fft_complex_high;
extern fftw_plan fft_forward;
extern fftw_plan fft_backwards;
extern fftw_plan fft_backwards_mid;
extern fftw_plan fft_backwards_high;


void fft_init();
void fft_clear();

void tlwe_sk_init_fft(tlwe_sk_fft *sk);
tlwe_sk_fft tlwe_keygen_fft(double param);
void tlwe_sk_clear_fft(tlwe_sk_fft sk);

void tgsw_sample_init_fft(tgsw_sample_fft *ct);
void tgsw_sample_clear_fft(tgsw_sample_fft ct);

void multiply_accumulate_poly_fft(uint64_t *out, uint64_t *a, fftw_complex *b);

tlwe_sample tlwe_encrypt_fft(tlwe_sk_fft sk, double param, int M, int* m);
void tlwe_encrypt_over_fft(tlwe_sample ct, tlwe_sk_fft sk, double param, int M, int* m);
tlwe_sample tlwe_encrypt_zero_fft(tlwe_sk_fft sk, double param);
void tlwe_encrypt_zero_over_fft(tlwe_sample ct, tlwe_sk_fft sk, double param);
tlwe_sample tlwe_encrypt_zero_gaussian_fft(tlwe_sk_fft sk, gaussian_param_t param);
void tlwe_encrypt_zero_gaussian_over_fft(tlwe_sample ct, tlwe_sk_fft sk, gaussian_param_t param);

int* tlwe_decrypt_fft(tlwe_sk_fft sk, int M, tlwe_sample ct);
void tlwe_decrypt_over_and_keep_fft(int *m, tlwe_sk_fft sk, int M, tlwe_sample ct, double *err);

lwe_sk tlwe_key_extract_fft(tlwe_sk_fft in);

pkc_fft sanitize_pkc_gen_fft(tlwe_sk_fft tsk, double parame);
void sanitize_pkc_clear_fft(pkc_fft PK);
void sanitize_pkc_enc_fft(tlwe_sample ct, pkc_fft PK, gaussian_param_t paramr, gaussian_param_t paramep, gaussian_param_t paramepp);

tgsw_sample_fft tgsw_encrypt_fft(tlwe_sk_fft sk, double param, int* m);
void tgsw_encrypt_over_fft(tgsw_sample_fft ct, tlwe_sk_fft sk, double param, int* m);

int* tgsw_decrypt_fft(tlwe_sk_fft sk, tgsw_sample_fft ct);
void tgsw_decrypt_over_and_keep_fft(int *m, tlwe_sk_fft sk, tgsw_sample_fft ct);

bsk_fft generate_bsk_fft(lwe_sk sk_in, tlwe_sk_fft sk_out, double param, int n);
void bsk_clear_fft(bsk_fft bsk, int n);

tlwe_sample randomized_external_product_fft(tlwe_sample tlwe, tgsw_sample_fft tgsw);
tlwe_sample blind_rotate_fft(lwe_sample lwe, int n, uint64_t* test_vector, bsk_fft bsk);
tlwe_sample cp_blind_rotate_fft(lwe_sample lwe, gaussian_param_t param, int n, uint64_t* testv, bsk_fft bsk);
//tlwe_sample cps_blind_rotate_online_fft(lwe_sample lwe, double param, double noiseparam, int M, int n, uint64_t* test_vector, bsk_fft bsk, tlwe_sk_fft sk);
void bootstrap_fft(lwe_sample ct, bsk_fft bsk, ksk ksk, int M_out, int n);
//void sanitize_s_fft(lwe_sample ct, bsk_fft bsk, ksk ksk, pks PK, double param, int M_out, int n);
void sanitize_c_fft(lwe_sample ct, bsk_fft bsk, ksk ksk, pkc_fft PK, gaussian_param_t paramy, gaussian_param_t paramr, gaussian_param_t paramep, gaussian_param_t paramepp, int M_out, int n);

/* DEBUG ONLY */
void randomized_accum_external_product_fft(tlwe_sample out, tlwe_sample tlwe, tgsw_sample_fft tgsw);
void accum_external_product_fft(tlwe_sample out, tlwe_sample tlwe, tgsw_sample_fft tgsw);

//rerand sanitize_pkenc_gen_fft(tlwe_sk_fft tsk, double param);
//tlwe_sample cp_blind_rotate_pkenc_fft(lwe_sample lwe, double param, int M, size_t n, uint64_t* test_vector, bsk_fft bsk, rerand PKenc);
//tlwe_sample cp_blind_rotate_pkenc_online_fft(lwe_sample lwe, double param, double noiseparam, int M, size_t n, uint64_t* test_vector, bsk_fft bsk, tlwe_sk_fft sk);
//void bootstrap_pkenc_fft(lwe_sample ct, bsk_fft bsk, ksk ksk, double param, int M_out, size_t n);
//void sanitize_pkenc_fft(lwe_sample ct, bsk_fft bsk, ksk ksk, rerand PKenc, double param, int M_out, size_t n);
#endif
