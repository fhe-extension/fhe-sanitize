//#pragma once
#include <stdio.h>
#include <stdlib.h>

#include "../fft.h"

#define _split (_ell * _logBg) / 3
#define _Ntries 1000

int main()
{
  fft_init();
  uint64_t *r = malloc(_N * sizeof(uint64_t));

  printf("Split : %d\n", _split);

  for (int i = 0; i < _N; ++i)
    r[i] = 0;
  r[0] = 3 << (64 - _logBg * _ell);
  r[1] = 4 << (64 - _logBg * _ell);
  uniform64_distribution_vector(r, 1);
  r[0] <<= 64 - _logBg * _ell;
  uint64_t c0 = r[0];
  printf("1st polynomial, coeff 0 = %llu\n", r[0]);

  fftw_complex *c1 = fftw_malloc(3 * (_N + 1) * sizeof(fftw_complex));

  for (int coeff = 0; coeff < _N; ++coeff){
    fft_real[coeff] = r[coeff] >> (64 - _split); // High order bits
  }
  memset(fft_real+_N, 0, _N * sizeof(double));
  fftw_execute(fft_forward);
  memcpy(c1, fft_complex, (_N+1) * sizeof(fftw_complex));
  for (int coeff = 0; coeff < _N; ++coeff){
    fft_real[coeff] = (r[coeff] << _split) >> (64 - _split); // Mid order bits
  }
  memset(fft_real+_N, 0, _N * sizeof(double));
  fftw_execute(fft_forward);
  memcpy(c1 + (_N+1), fft_complex, (_N+1) * sizeof(fftw_complex));
  for (int coeff = 0; coeff < _N; ++coeff){
    fft_real[coeff] = (r[coeff] << (2 * _split)) >> (64 - _split); // Low order bits
  }
  memset(fft_real+_N, 0, _N * sizeof(double));
  fftw_execute(fft_forward);
  memcpy(c1 + 2 * (_N+1), fft_complex, (_N+1) * sizeof(fftw_complex));

  for (int i = 0; i < _N; ++i)
    r[i] = 0;
  r[0] = 1;
  r[1] = 2;
  uniform64_distribution_vector(r, 1);
  printf("2nd polynomial, coeff 0 = %llu\n", r[0]);


  printf("Expected coeff 0 : %llu\n", c0 * r[0]);
  uint64_t out[_N];
  for (int i = 0; i < _N; ++i)
    out[i] = 0;

  multiply_accumulate_poly_fft(out, r, c1);
  for (int i = 0; i < 10; ++i)
    printf("%llu\n", out[i]);

  /* FULL SPACE PRODUCT */
  int cpt = 0;
  double var = 0;
  printf("FULL SPACE PRODUCT TEST\n");
  for (int i = 0; i < _Ntries; ++i){
    uniform64_distribution_vector(r, _N);
    for (int i = 0; i < _N; ++i)
      r[i] <<= 64 - _logBg * _ell;

    for (int coeff = 0; coeff < _N; ++coeff){
      fft_real[coeff] = r[coeff] >> (64 - _split); // High order bits
    }
    memset(fft_real+_N, 0, _N * sizeof(double));
    fftw_execute(fft_forward);
    memcpy(c1, fft_complex, (_N+1) * sizeof(fftw_complex));
    for (int coeff = 0; coeff < _N; ++coeff){
      fft_real[coeff] = (r[coeff] << _split) >> (64 - _split); // Mid order bits
    }
    memset(fft_real+_N, 0, _N * sizeof(double));
    fftw_execute(fft_forward);
    memcpy(c1 + (_N+1), fft_complex, (_N+1) * sizeof(fftw_complex));
    for (int coeff = 0; coeff < _N; ++coeff){
      fft_real[coeff] = (r[coeff] << (2 * _split)) >> (64 - _split); // Low order bits
    }
    memset(fft_real+_N, 0, _N * sizeof(double));
    fftw_execute(fft_forward);
    memcpy(c1 + 2 * (_N+1), fft_complex, (_N+1) * sizeof(fftw_complex));

    uint64_t r2[_N];
    uniform64_distribution_vector(r2, _N);

    uint64_t out2[_N];
    for (int i = 0; i < _N; ++i){
      r2[i] &= (1L << _logBg * _ell) - 1;
      out[i] = 0;
      out2[i] = 0;
    }

    multiply_accumulate_poly_fft(out, r2, c1);

    multiply_accumulate_poly(out2, r2, r);
    for (int i = 0; i < _N; ++i){
      if (out[i] != out2[i]){
        cpt++;
        printf("MISMATCH %d : %llu == %llu --- ", i, out[i], out2[i]);
        printf("Delta : %lld\n", ((int64_t) out[i] - (int64_t) out2[i]));
      }
      var += ((int64_t) out[i] - (int64_t) out2[i]) * ((int64_t) out[i] - (int64_t) out2[i]);
    }
  }

  printf("Mismatchs : %d; Variance : %lf\n", cpt, var / _N);

  fft_clear();
  free(r);
  free(c1);

#ifdef _WIN32
  system("pause"); // Pauses to actually see the output in windows
#endif

  return 0;
}


