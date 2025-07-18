#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

#include "../fft.h"
//#include "../tlwe.h"
#include "../clock.h"

#define TRIES 1e4


#define N_TRIES 1000

#define TEST_ENC_DEC 0

#define _two64_double powl(2.0L,64.0L)
#define _PI acosl(-1.0L)

int main(int argc, char const *argv[])
{
  srand(time(NULL));
  fft_init();

  double q = powl(2.0L,_ell*_logBg);
  double var_cdot = powl(2.0L,-41.7L);
  gaussian_param_t param_boxdot = gaussian(sqrt(2.0L*_PI*var_cdot)*q);
  double stddev_bk = powl(2.0L,-42.0L);
  precompute_random(param_boxdot);

  double accum_rnd = 0, accum_enc = 0, accum_dec = 0;
  double param = (uint64_t) 1 << (64 - _logBg - 4);

  tlwe_sk_fft sk = tlwe_keygen_fft(1.0L);
  lwe_sk ext_key = tlwe_key_extract_fft(sk);
  lwe_sample ext_ct;
  lwe_sample_init(&ext_ct,_N*_k);


  tgsw_sample_fft ct;
  tgsw_sample_init_fft(&ct);
  tlwe_sample ctt;
  tlwe_sample_init(&ctt);

  int* m= (int*)malloc(_N * sizeof(int));
  int* one = (int*) malloc(_N * sizeof(int));

  for(int j = 0; j < _N; ++j){
    m[j] = rand() % _Bg;
    one[j] = 0;
  }
  one[0] = 1;

  double var = 0.0;
  for(int try = 0; try < N_TRIES; ++try) {

    tgsw_encrypt_over_fft(ct, sk, stddev_bk*_two64_double, one);
    tlwe_encrypt_over_fft(ctt, sk, stddev_bk*_two64_double, _Bg, m);

    tlwe_sample out;
    tlwe_sample_init(&out);
    for (int i = 0; i < (_k + 1) * _N; ++i)
      out[i] = 0;

    start_chrono();
    accum_external_product_fft(out, ctt, ct);
    accum_enc += stop_chrono();
    accum_dec += stop_chrono();


    tlwe_extract_over_and_keep(ext_ct, out);
    int *mm = tlwe_decrypt_fft(sk, _Bg, out);



    // compare the polynomials
    for (int i = 0; i < _N; ++i) {
      //printf("%d == %d\n", m[i] % _Bg, mm[i] % _Bg);
      assert(m[i] % _Bg == mm[i] % _Bg);
    }
    free(mm);

    // measuring error variance
    double err;
    lwe_decrypt_and_keep(ext_key, _Bg, ext_ct, _k * _N, &err);
    var += err * err;

  }

  printf("Time for generating a random polynomial: %f microseconds\n", (double) accum_rnd / N_TRIES);
  printf("Time for encryption:                     %f microseconds\n", (double) accum_enc / N_TRIES);
  printf("Time for decryption:                     %f microseconds\n", (double) accum_dec / N_TRIES);
  printf("Measured variance for multiplying by TGSW(1) : %le\n", var / N_TRIES);

  var = 0.0;
  one[0] = 0;
  for(int try = 0; try < N_TRIES; ++try) {

    tgsw_encrypt_over_fft(ct, sk, stddev_bk*_two64_double, one);
    tlwe_encrypt_over_fft(ctt, sk, stddev_bk*_two64_double, _Bg, m);

    tlwe_sample out;
    tlwe_sample_init(&out);
    for (int i = 0; i < (_k + 1) * _N; ++i)
      out[i] = 0;

    start_chrono();
    accum_external_product_fft(out, ctt, ct);
    accum_enc += stop_chrono();
    accum_dec += stop_chrono();


    tlwe_extract_over_and_keep(ext_ct, out);
    int *mm = tlwe_decrypt_fft(sk, _Bg, out);



    // compare the polynomials
    for (int i = 0; i < _N; ++i) {
      //printf("0 == %d\n", mm[i] % _Bg);
      assert(0 == mm[i] % _Bg);
    }
    free(mm);


    // measuring error variance
    double err;
    lwe_decrypt_and_keep(ext_key, _Bg, ext_ct, _k * _N, &err);
    var += err * err;

  }
  printf("Measured variance for multiplying by TGSW(0) : %le\n", var / N_TRIES);

  // cleanup
  free(m);
  free(one);
  tlwe_sk_clear_fft(sk);
  lwe_sk_clear(ext_key);
  lwe_sample_clear(ext_ct);
  tlwe_sample_clear(ctt);
  tgsw_sample_clear_fft(ct);

#ifdef _WIN32
  system("pause"); // Pauses to actually see the output in windows
#endif

  return 0;
}
