#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>

#include "../clock.h"
#include "../fft.h"

#define _two64_double powl(2.0L,64.0L)
#define _PI acosl(-1.0L)

int main(int argc, char const *argv[])
{
  srand(time(NULL));
  fft_init();

  const int MSG_SPACE = 2;
  const int _n = 538;

  long double stddev_bk = powl(2.0L,-33.8L);
  long double stddev_ks = powl(2.0L,-13.6L);
  long double var_cdot = powl(2.0L,-39.4L);
  long double var_r = powl(2.0L,-53.2L);
  long double var_e2 = powl(2.0L,-38.1L);

  //long double invq = powl(2.0L,64.0L-_ell*_logBg);
  long double q = powl(2.0L,_ell*_logBg);

  gaussian_param_t param_boxdot = gaussian(sqrtl(2.0L*_PI*var_cdot)*q);
  gaussian_param_t param_r = gaussian(sqrtl(2*_PI*var_r)*q);
  gaussian_param_t param_e2 = gaussian(sqrtl(2*_PI*var_e2)*q);
  precompute_random(param_boxdot);

  lwe_sk ns = lwe_keygen(_n);
  tlwe_sk_fft S = tlwe_keygen_fft(stddev_bk*q);
  lwe_sk Ns = tlwe_key_extract_fft(S);

  ksk ksk = generate_ksk(Ns, ns, stddev_ks*_two64_double, _k * _N, _n);
  bsk_fft bsk = generate_bsk_fft(ns, S, stddev_bk*_two64_double, _n);

  printf("Precompute PkEnc\n");

  pkc_fft PK = sanitize_pkc_gen_fft(S, sqrtl(2*_PI)*stddev_bk*q);

  int m,m2;
  m = rand()% 2;

  lwe_sample ct = lwe_encrypt(ns, MSG_SPACE, powl(2.0L,50.0L), m, _n);

  printf("GO !\n");

  start_chrono();

  sanitize_c_fft(ct, bsk, ksk, PK, param_boxdot, param_r, param_r, param_e2, MSG_SPACE,_n);
  //bootstrap_fft(ct, bsk, ksk, MSG_SPACE,_n);


  long double temps = stop_chrono();

  printf("STOP !\n");

  m2 = lwe_decrypt_and_keep(ns, MSG_SPACE, ct,_n);

  printf("%u == %u\n", m, m2);

  printf("Time for sanitization: %lf microseconds\n", (double) temps);

  lwe_sample_clear(ct);
  lwe_sk_clear(Ns);
  lwe_sk_clear(ns);
  ksk_clear(ksk, _k * _N, _n);
  sanitize_pkc_clear_fft(PK);
  bsk_clear_fft(bsk, _n);
  tlwe_sk_clear_fft(S);
  clear_random();

  fft_clear();

#ifdef _WIN32
  system("pause"); // Pauses to actually see the output in windows
#endif

  return 0;
}
