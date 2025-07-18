//#pragma once

#include "tlwe.h"
#include <string.h>


#define _two64_double powl(2.0L,64.0L)

void tlwe_sample_init(tlwe_sample *ct)
{
  *ct = (tlwe_sample) malloc((_k+1) * _N * sizeof(uint64_t));
}

void tlwe_sample_clear(tlwe_sample ct)
{
  free(ct);
}

void tlwe_sk_init(tlwe_sk *sk)
{
  *sk = (tlwe_sk) malloc (_k * _N * sizeof(uint64_t));
}

void tlwe_sk_clear(tlwe_sk sk)
{
  free(sk);
}

void tlwe_sample_zero(tlwe_sample ct)
{
  int max = (_k+1)*_N;
  for(int i=0; i < max; ++i) {
    ct[i] = 0;
  }
}

tlwe_sk tlwe_keygen()
{
  tlwe_sk tsk;
  tlwe_sk_init(&tsk);
  random_binary_vector(tsk, _k*_N);
  return tsk;
}

// out += a
void tlwe_add_to(tlwe_sample out, tlwe_sample a)
{
  for (int i = 0; i < (_k + 1) * _N; ++i)
    out[i] += a[i];
}

// out += a * b // Inputs are degree N poly modulo X^N+1
void multiply_accumulate_poly(uint64_t *out, uint64_t *a, uint64_t *b)
{
  int coeff1, coeff2;
  for (coeff1 = 0; coeff1 < _N; ++coeff1)
    {
      for (coeff2 = 0; coeff2 <= coeff1; coeff2++)
        out[coeff1] += a[coeff2]*b[coeff1-coeff2];
      for (; coeff2 < _N; coeff2++)
        out[coeff1] -= a[coeff2]*b[_N+coeff1-coeff2];
    }
}

void tlwe_encrypt_over(tlwe_sample ct, tlwe_sk sk, double param, int M, int *m)
{
  //tlwe_sample_zero(ct, k, N);

  int max=_k*_N;

  //noise polynomial e
  double *e = (double*) malloc(_N * sizeof(double));

  //b de taille N b=e gaussien
  noise_vector(e,param,_N);

  int poly, coeff1, coeff2;

  for (coeff1 = 0; coeff1 < _N; ++coeff1)
    ct[max+coeff1] = ((uint64_t) (e[coeff1] + (_two64_double* m[coeff1])/M)) & _mask;
  free(e);

  uniform64_distribution_vector(ct, max);
  for (poly = 0; poly < _k; ++poly)
    for (coeff1 = 0; coeff1 < _N; ++coeff1)
      ct[poly*_N+coeff1] &= _mask;


  for (poly = 0; poly < _k; ++poly)
    //multiply_accumulate_poly(ct+max,ct+_N*poly,sk+_N*poly);
    for (coeff1 = 0; coeff1 < _N; ++coeff1)
      {
        for (coeff2 = 0; coeff2 <= coeff1; coeff2++)
          ct[max+coeff1] += ct[poly*_N+coeff2]*sk[poly*_N+coeff1-coeff2];
        for (; coeff2 < _N; coeff2++)
          ct[max+coeff1] -= ct[poly*_N+coeff2]*sk[poly*_N+_N+coeff1-coeff2];
      }
}
tlwe_sample tlwe_encrypt(tlwe_sk sk, double param, int M, int *m)
{
  tlwe_sample ct;
  tlwe_sample_init(&ct);

  tlwe_encrypt_over(ct, sk, param, M, m);

  return ct;
}

void tlwe_encrypt_zero_over(tlwe_sample ct, tlwe_sk sk, double param)
{
  int max=_k*_N;

  //noise polynomial e
  double *e = (double*) malloc(_N * sizeof(double));

  noise_vector(e,param,_N);

  int poly, coeff1, coeff2;

  for (coeff1 = 0; coeff1 < _N; ++coeff1)
    ct[max+coeff1] = ((uint64_t) e[coeff1]) & _mask;
  free(e);

  uniform64_distribution_vector(ct, max);
  for (poly = 0; poly < _k; ++poly)
    for (coeff1 = 0; coeff1 < _N; ++coeff1)
      ct[poly*_N+coeff1] &= _mask;

  for (poly = 0; poly < _k; ++poly)
    for (coeff1 = 0; coeff1 < _N; ++coeff1)
      {
        for (coeff2 = 0; coeff2 <= coeff1; coeff2++)
          ct[max+coeff1] += ct[poly*_N+coeff2]*sk[poly*_N+coeff1-coeff2];
        for (; coeff2 < _N; coeff2++)
          ct[max+coeff1] -= ct[poly*_N+coeff2]*sk[poly*_N+_N+coeff1-coeff2];
      }
}
tlwe_sample tlwe_encrypt_zero(tlwe_sk sk, double param)
{
  tlwe_sample ct;
  tlwe_sample_init(&ct);

  tlwe_encrypt_zero_over(ct, sk, param);

  return ct;
}

void tlwe_encrypt_zero_gaussian_over(tlwe_sample ct, tlwe_sk sk, gaussian_param_t param)
{
  int max=_k*_N;

  int poly, coeff1, coeff2;

  gaussian_overZ_vector(ct+max,param,_N);
  for (coeff1 = 0; coeff1 < _N; ++coeff1)
    ct[max+coeff1] <<= (64-_logBg*_ell);

  uniform64_distribution_vector(ct, max);
  for (poly = 0; poly < _k; ++poly)
    for (coeff1 = 0; coeff1 < _N; ++coeff1)
      ct[poly*_N+coeff1] &= _mask;

  for (poly = 0; poly < _k; ++poly)
    for (coeff1 = 0; coeff1 < _N; ++coeff1)
      {
        for (coeff2 = 0; coeff2 <= coeff1; coeff2++)
          ct[max+coeff1] += ct[poly*_N+coeff2]*sk[poly*_N+coeff1-coeff2];
        for (; coeff2 < _N; coeff2++)
          ct[max+coeff1] -= ct[poly*_N+coeff2]*sk[poly*_N+_N+coeff1-coeff2];
      }
}
tlwe_sample tlwe_encrypt_zero_gaussian(tlwe_sk sk, gaussian_param_t param)
{
  tlwe_sample ct;
  tlwe_sample_init(&ct);

  tlwe_encrypt_zero_gaussian_over(ct, sk, param);

  return ct;
}

int* tlwe_decrypt(tlwe_sk sk, int M, tlwe_sample ct)
{
  int max=_k*_N;

  int *m = (int *) malloc(_N * sizeof(int));

  int poly, coeff1, coeff2;
  for (poly = 0; poly < _k; ++poly)
    for (coeff1 = 0; coeff1 < _N; ++coeff1)
      {
        for (coeff2 = 0; coeff2 <= coeff1; coeff2++)
          ct[max+coeff1] -= ct[poly*_N+coeff2]*sk[poly*_N+coeff1-coeff2];
        for (; coeff2 < _N; coeff2++)
          ct[max+coeff1] += ct[poly*_N+coeff2]*sk[poly*_N+_N+coeff1-coeff2];
      }
  for (coeff1 = 0; coeff1 < _N; ++coeff1)
    m[coeff1]=(int)(((double) M * ct[max+coeff1])/_two64_double + 0.5);
  tlwe_sample_clear(ct);
  return m;
}

void tlwe_decrypt_over_and_keep(int* m, tlwe_sk sk, int M, tlwe_sample ct)
{
  int max=_k*_N;

  uint64_t *b = (uint64_t *) malloc(_N * sizeof(uint64_t));
  memcpy(b, ct+max, _N * sizeof(uint64_t));

  int poly, coeff1, coeff2;
  for (poly = 0; poly < _k; ++poly)
    for (coeff1 = 0; coeff1 < _N; ++coeff1)
      {
        for (coeff2 = 0; coeff2 <= coeff1; coeff2++)
          b[coeff1] -= ct[poly*_N+coeff2]*sk[poly*_N+coeff1-coeff2];
        for (; coeff2 < _N; coeff2++)
          b[coeff1] += ct[poly*_N+coeff2]*sk[poly*_N+_N+coeff1-coeff2];
      }
  for (coeff1 = 0; coeff1 < _N; ++coeff1)
    m[coeff1]=(int)(((double) M * b[coeff1])/_two64_double + 0.5);
  free(b);
}


lwe_sample tlwe_extract(tlwe_sample in)
{
  lwe_sample out;
  lwe_sample_init(&out, _k * _N);
  out[_k*_N] = in[_k*_N];

  int i,j;
  int a_pos = 0;
  for(i = 0; i < _k; ++i) {
    out[a_pos++] = in[_N*i];
    for(j = _N - 1; j > 0; --j) {
      out[a_pos++] = (-in[_N*i+j]);
    }
  }
  tlwe_sample_clear(in);
  return out;
}

void tlwe_extract_over_and_keep(lwe_sample out, tlwe_sample in)
{
  int a_pos = 0;
  //nieme position pour le coeff n+1 --> coeff constant de in[]
  out[_k*_N] = in[_k*_N];

  //taille de la clé =k*N=n-1
  //tlwe=(a',b')--> lwe=(a,b) witj b=b'_0 constant term of b'
  // coeff(a)_k = coeff(a'_k)/X (décallé )
  for(int i = 0; i < _k; ++i) {
    out[a_pos++] = in[_N*i];

    // reverse and change the sign of the others
    for(int j = _N - 1; j > 0; --j) {
      out[a_pos++] = (-in[_N*i+j]);
    }
  }
}

//n=kN
lwe_sk tlwe_key_extract(tlwe_sk in)
{
  int sk_pos = 0;
  lwe_sk sk_out = (uint64_t *) malloc(_k * _N * sizeof(uint64_t));
  // k*N=n-1
  // concatenate the coefficients of the secret key polynomials
  //of size k*N
  for(int i = 0; i < _k; ++i) {
    for(int j = 0; j < _N; ++j) {
      sk_out[sk_pos++] = in[i*_N+j];
    }
  }
  return sk_out;
}


//B and ell hard coded
//return a unint64* table of size (N*ell)*(k+1) = e_11, ..., ek+1,ell in R_q
//s.t.
void decompose_poly(uint64_t *out, uint64_t *in)
{
  int coeff, i;
  for (coeff = 0; coeff < _N; ++coeff)
    for (i=0; i < _ell; ++i)
      out[i*_N+coeff] = (in[coeff] << i * _logBg) >> (64 - _logBg);
}
void decompose_tlwe(uint64_t *out, tlwe_sample in)
{
  for (int i = 0; i < _k+1; ++i)
    decompose_poly(out+(i*_ell*_N), in + i*_N);
}
