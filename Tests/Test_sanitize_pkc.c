#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <fcntl.h> // for open
#include <unistd.h> // for close

#include "../clock.h"
#include "../fft.h"
//#include "../lwe.h"
//#include "../random.h"
//#include "Test_parameters.h"

#define N_TRIES 1000

#define TLWE_MSG_SPACE 4

#define _two64_double powl(2.0L,64.0L)
#define _PI acosl(-1.0L)


int main(int argc, char const *argv[])
{
  double accum_sanitize = 0;
  srand(time(NULL));

  double q = powl(2.0L,_ell*_logBg);

  double stddev_bk = powl(2.0L,-33.8L);
  //double var_r = 4*stddev_bk*stddev_bk;
  double var_e2 = 4*stddev_bk*stddev_bk*2*_PI*q*q*stddev_bk*stddev_bk*_N*2;
  double var_r = powl(2.0L,-62.6L);
  //double var_r = powl(2.0L,-53.2L);
  //double var_e2 = powl(2.0L,-51L);

  gaussian_param_t param_r = gaussian(sqrtl(2*_PI*var_r)*q);
  gaussian_param_t param_e2 = gaussian(sqrtl(2*_PI*var_e2)*q);


  fft_init();

  tlwe_sk_fft tsk = tlwe_keygen_fft(stddev_bk*q);
  tlwe_sample ct;
  tlwe_sample_init(&ct);

  pkc_fft PK = sanitize_pkc_gen_fft(tsk, sqrtl(2*_PI)*stddev_bk*q);

  double *err = malloc(_N * sizeof(double));
  double accum_var = 0.;

  printf("PK generated\n");

#ifdef _WIN32
  system("pause"); // Pauses to actually see the output in windows
#endif

  int* mu= (int*)malloc(_N * sizeof(int));
  int* mmu= (int*)malloc(_N * sizeof(int));
  for (int try = 0; try < N_TRIES; ++try){
    size_t j;
    for(j = 0; j < _N; ++j)
      mu[j]=  rand()% TLWE_MSG_SPACE;

    tlwe_encrypt_over_fft(ct, tsk, 0, TLWE_MSG_SPACE, mu);

    //tlwe_sample_zero(ct);

    //printf("Sanitizing\n");

    start_chrono();

    sanitize_pkc_enc_fft(ct, PK, param_r, param_r, param_e2);

    accum_sanitize += stop_chrono();
    //printf("Sanitized\n");
    tlwe_decrypt_over_and_keep_fft(mmu, tsk, TLWE_MSG_SPACE, ct, err);

    for (int coeff = 0; coeff < _N; ++coeff)
      accum_var += err[coeff] * err[coeff] * powl(10.0L,12.0L);

    for (j = 0; j < _N; ++j)
      if (mu[j] % TLWE_MSG_SPACE != mmu[j] % TLWE_MSG_SPACE)
        printf("%u == %u\n",mu[j],mmu[j]);
    //printf("\n");
  }
  printf("Time for adding a public key encryption of 0 to a TLWE sample: %f microseconds\n", (double) accum_sanitize/N_TRIES);
  printf("Variance of the error : %f\n", (double)accum_var/(_N*N_TRIES));

  free(err);
  free(mu);
  free(mmu);
  //sanitize_pkc_clear_fft(PK);
  tlwe_sk_clear_fft(tsk);
  tlwe_sample_clear(ct);

  fft_clear();

#ifdef _WIN32
  system("pause"); // Pauses to actually see the output in windows
#endif

  return 0;
}
