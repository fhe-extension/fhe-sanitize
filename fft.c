
#include "fft.h"

#define _two64_double powl(2.0L,64.0L)

#define _split ((_ell * _logBg) / 3)

double *fft_real;
fftw_complex *fft_complex;
fftw_complex *fft_complex_mid;
fftw_complex *fft_complex_high;
fftw_plan fft_forward;
fftw_plan fft_backwards;
fftw_plan fft_backwards_mid;
fftw_plan fft_backwards_high;



void fft_init()
{
  fft_real = (double *) fftw_malloc(sizeof(double) * 2 * _N);

  fft_complex = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (_N+1));
  fft_complex_mid = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (_N+1));
  fft_complex_high = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (_N+1));

  fftw_import_wisdom_from_filename("wisdom");
  fft_forward = fftw_plan_dft_r2c_1d(2 * _N, fft_real, fft_complex, FFTW_EXHAUSTIVE | FFTW_DESTROY_INPUT);
  fft_backwards = fftw_plan_dft_c2r_1d(2 * _N, fft_complex, fft_real, FFTW_EXHAUSTIVE | FFTW_DESTROY_INPUT);
  fft_backwards_mid = fftw_plan_dft_c2r_1d(2 * _N, fft_complex_mid, fft_real, FFTW_EXHAUSTIVE | FFTW_DESTROY_INPUT);
  fft_backwards_high = fftw_plan_dft_c2r_1d(2 * _N, fft_complex_high, fft_real, FFTW_EXHAUSTIVE | FFTW_DESTROY_INPUT);
  fftw_export_wisdom_to_filename("wisdom");
}
void fft_clear()
{
  fftw_destroy_plan(fft_forward); fftw_destroy_plan(fft_backwards);
  fftw_free(fft_real); fftw_free(fft_complex);
}

void tlwe_sk_init_fft(tlwe_sk_fft *sk)
{
  *sk = (tlwe_sk_fft) fftw_malloc(3 * _k * (_N+1) * sizeof(fftw_complex));
}
tlwe_sk_fft tlwe_keygen_fft(double param)
{
  double *key = (double *) malloc(_N * sizeof(double));
  uint64_t *nofft = (uint64_t *) malloc(_N * sizeof(uint64_t));
  tlwe_sk_fft tsk;
  tlwe_sk_init_fft(&tsk);
  for (int i = 0; i < _k; ++i){
    noise_vector(key, param, _N);
    double tmp;
    for (int coeff = 0; coeff < _N; ++coeff){
    tmp= key[coeff] < 0 ? (key[coeff] + _two64_double) : key[coeff];
    nofft[coeff] = fmodl(round(tmp),_two64_double);
    }
    for (int coeff = 0; coeff < _N; ++coeff)
      fft_real[coeff] = nofft[coeff] >> (64 - _split); // High order bits
    memset(fft_real+_N, 0, _N * sizeof(double));
    fftw_execute(fft_forward);
    memcpy(tsk + 3 * i * (_N+1), fft_complex, (_N+1) * sizeof(fftw_complex));
    for (int coeff = 0; coeff < _N; ++coeff)
      fft_real[coeff] = (nofft[coeff] << _split) >> (64 - _split); // Mid order bits
    memset(fft_real+_N, 0, _N * sizeof(double));
    fftw_execute(fft_forward);
    memcpy(tsk + (3 * i + 1) * (_N+1), fft_complex, (_N+1) * sizeof(fftw_complex));
    for (int coeff = 0; coeff < _N; ++coeff)
      fft_real[coeff] = (nofft[coeff] << (2 * _split)) >> (64 - _split); // Low order bits
    memset(fft_real+_N, 0, _N * sizeof(double));
    fftw_execute(fft_forward);
    memcpy(tsk + (3 * i + 2) * (_N+1), fft_complex, (_N+1) * sizeof(fftw_complex));
  }
  free(key);
  free(nofft);
  return tsk;
}
void tlwe_sk_clear_fft(tlwe_sk_fft sk)
{
  fftw_free(sk);
}

void tgsw_sample_init_fft(tgsw_sample_fft *ct)
{
  *ct = (tgsw_sample_fft) fftw_malloc(3 * (_k+1) * (_k+1) * (_N+1) * _ell * sizeof(fftw_complex));
}
void tgsw_sample_clear_fft(tgsw_sample_fft ct)
{
  fftw_free(ct);
}

void multiply_accumulate_poly_fft(uint64_t *out, uint64_t *a, fftw_complex *b)
{
  for (int i = 0; i < _N; ++i)
    fft_real[i] = (a[i] << (64 - _split)) >> (64 - _split); // Low order bits
  memset(fft_real+_N, 0, _N * sizeof(double));
  fftw_execute(fft_forward);
  for (int i = 0; i <= _N; ++i){
    fft_complex_high[i] = fft_complex[i] * b[i]; // * high order bits from bsk
    fft_complex_mid[i] = fft_complex[i] * b[i + (_N + 1)]; // * mid order bits from bsk
    fft_complex[i] *= b[i + 2 * (_N + 1)]; // * low order bits from bsk
  }
  fftw_execute(fft_backwards);
  for (int i = 0; i < _N; ++i){
    double tmpf = (fft_real[i] + _N) / (2 * _N);
    //assert(0 <= tmpf && tmpf < _two64_double);
    uint64_t tmp = tmpf;
    tmp <<= 64 - 3 * _split; // Low * low
    out[i] += tmp;
    tmpf = (fft_real[_N + i] + _N) / (2 * _N);
    //assert(0 <= tmpf && tmpf < _two64_double);
    tmp = tmpf;
    tmp <<= 64 - 3 * _split; // Low * low
    out[i] -= tmp;
    //assert(0 <= out[i] && out[i] < _two64_double);
  }
  for (int i = 0; i < _N; ++i)
    fft_real[i] = (a[i] << (64 - 2 * _split)) >> (64 - _split); // Mid order bits
  memset(fft_real+_N, 0, _N * sizeof(double));
  fftw_execute(fft_forward);
  for (int i = 0; i <= _N; ++i){
    fft_complex_high[i] += fft_complex[i] * b[i + (_N + 1)]; // * mid order bits from bsk
    fft_complex_mid[i] += fft_complex[i] * b[i + 2 * (_N + 1)]; // * low order bits from bsk
  }
  fftw_execute(fft_backwards_mid);
  for (int i = 0; i < _N; ++i){
    double tmpf = (fft_real[i] + _N) / (2 * _N);
    //assert(0 <= tmpf && tmpf < _two64_double);
    uint64_t tmp = tmpf;
    tmp <<= 64 - 2 * _split; // Low * mid + Mid * low
    out[i] += tmp;
    tmpf = (fft_real[_N + i] + _N) / (2 * _N);
    //assert(0 <= tmpf && tmpf < _two64_double);
    tmp = tmpf;
    tmp <<= 64 - 2 * _split; // Low * mid + Mid * low
    out[i] -= tmp;
  }
  for (int i = 0; i < _N; ++i)
    fft_real[i] = (a[i] << (64 - 3 * _split)) >> (64 - _split); // High order bits
  memset(fft_real+_N, 0, _N * sizeof(double));
  fftw_execute(fft_forward);
  for (int i = 0; i <= _N; ++i){
    fft_complex_high[i] += fft_complex[i] * b[i + 2 * (_N + 1)]; // * low order bits from bsk
  }
  fftw_execute(fft_backwards_high);
  for (int i = 0; i < _N; ++i){
    double tmpf = (fft_real[i] + _N) / (2 * _N);
    //assert(0 <= tmpf && tmpf < _two64_double);
    uint64_t tmp = tmpf;
    tmp <<= 64 - _split; // Low * high + Mid * mid + High * low
    out[i] += tmp;
    tmpf = (fft_real[_N + i] + _N) / (2 * _N);
    //assert(0 <= tmpf && tmpf < _two64_double);
    tmp = tmpf;
    tmp <<= 64 - _split; // Low * high + Mid * mid + High * low
    out[i] -= tmp;
  }
}

void tlwe_dotproduct_fft(tlwe_sample ct, tlwe_sk_fft sk){
  int max=_k*_N;
  for (int i = 0; i < _k; ++i)
    multiply_accumulate_poly_fft(ct + max, ct + i, sk + i);
}

void tlwe_encrypt_over_fft(tlwe_sample ct, tlwe_sk_fft sk, double param, int M, int* m)
{
  int max=_k*_N;

  //noise polynomial e
  double *e = (double*) malloc(_N * sizeof(double));

  noise_vector(e,param,_N);

  for (int i = 0; i < _N; ++i){
    double tmp = e[i] + (_two64_double * m[i]) / M;
    tmp = tmp < 0 ? tmp + _two64_double : tmp;
    ct[max+i] = round(tmp);
    ct[max+i] &= _mask;
  }
  free(e);

  uniform64_distribution_vector(ct, max);
  for (int i = 0; i < max; ++i)
    ct[i] &= _mask;

  tlwe_dotproduct_fft(ct, sk);

}
tlwe_sample tlwe_encrypt_fft(tlwe_sk_fft sk, double param, int M, int* m)
{
  tlwe_sample ct;
  tlwe_sample_init(&ct);
  tlwe_encrypt_over_fft(ct, sk, param, M, m);
  return ct;
}

void tlwe_encrypt_zero_over_fft(tlwe_sample ct, tlwe_sk_fft sk, double param)
{
  int max=_k*_N;

  //noise polynomial e
  double *e = (double*) malloc(_N * sizeof(double));

  noise_vector(e,param,_N);

  for (int i = 0; i < _N; ++i){
    double tmp = e[i] < 0 ? e[i] + _two64_double : e[i];
    //printf("%llu\n",tmp);
    ct[max+i] =  fmodl(round(tmp),_two64_double); 
    ct[max+i] &= _mask;
  }
  free(e);

  uniform64_distribution_vector(ct, max);
  for (int i = 0; i < max; ++i)
    ct[i] &= _mask;

  tlwe_dotproduct_fft(ct, sk);
}
tlwe_sample tlwe_encrypt_zero_fft(tlwe_sk_fft sk, double param)
{
  tlwe_sample ct;
  tlwe_sample_init(&ct);
  tlwe_encrypt_zero_over_fft(ct, sk, param);
  return ct;
}

void tlwe_encrypt_zero_gaussian_over_fft(tlwe_sample ct, tlwe_sk_fft sk, gaussian_param_t param)
{
  int max=_k*_N;

  gaussian_overZ_vector(ct+max,param,_N);
  for (int i = 0; i < _N; ++i)
    ct[max+i] <<= (64-_logBg*_ell);

  uniform64_distribution_vector(ct, max);
  for (int i = 0; i < max; ++i)
    ct[i] &= _mask;

  tlwe_dotproduct_fft(ct, sk);
}
tlwe_sample tlwe_encrypt_zero_gaussian_fft(tlwe_sk_fft sk, gaussian_param_t param)
{
  tlwe_sample ct;
  tlwe_sample_init(&ct);
  tlwe_encrypt_zero_gaussian_over_fft(ct, sk, param);
  return ct;
}

int* tlwe_decrypt_fft(tlwe_sk_fft sk, int M, tlwe_sample ct)
{
  int *m = (int *) malloc(_N * sizeof(int));
  double err[_N];
  tlwe_decrypt_over_and_keep_fft(m, sk, M, ct, err);
  tlwe_sample_clear(ct);
  return m;
}
void tlwe_decrypt_over_and_keep_fft(int *m, tlwe_sk_fft sk, int M, tlwe_sample ct, double *err)
{
  int max=_k*_N;

  uint64_t *b = (uint64_t *) malloc(_N * sizeof(uint64_t));
  for (int i = 0; i < _N; ++i)
    b[i] = -ct[max + i];

  for (int i = 0; i < _k; ++i)
    multiply_accumulate_poly_fft(b, ct + i, sk + i);
  for (int i = 0; i < _N; ++i){
    m[i]=(int)(((double) M * (-b[i]))/_two64_double + 0.5);
    err[i] = (((double) (-b[i])) - _two64_double * (double) m[i] / M) / _two64_double;
    m[i] = m[i] % M;
  }
  free(b);
}

lwe_sk tlwe_key_extract_fft(tlwe_sk_fft in)
{
  tlwe_sk nofft;
  tlwe_sk_init(&nofft);
  for (int i = 0; i < _k; ++i){
    memcpy(fft_complex, in + 3 * i * (_N + 1), (_N + 1) * sizeof(fftw_complex));
    fftw_execute(fft_backwards);
    for (int j = 0; j < _N; ++j){
      nofft[j + i * _N] = (fft_real[j] - fft_real[_N + j] + _N) / (2 * _N);
      nofft[j + i * _N] <<= (64 - _split);
    }
    memcpy(fft_complex, in + (3 * i + 1) * (_N + 1), (_N + 1) * sizeof(fftw_complex));
    fftw_execute(fft_backwards);
    for (int j = 0; j < _N; ++j){
      uint64_t tmp = (fft_real[j] - fft_real[_N + j] + _N) / (2 * _N);
      nofft[j + i * _N] += tmp << (64 - 2 * _split);
    }
    memcpy(fft_complex, in + (3 * i + 2) * (_N + 1), (_N + 1) * sizeof(fftw_complex));
    fftw_execute(fft_backwards);
    for (int j = 0; j < _N; ++j){
      uint64_t tmp = (fft_real[j] - fft_real[_N + j] + _N) / (2 * _N);
      nofft[j + i * _N] += tmp << (64 - 3 * _split);
    }
  }
  lwe_sk out = tlwe_key_extract(nofft);
  tlwe_sk_clear(nofft);
  return out;
}

pkc_fft sanitize_pkc_gen_fft(tlwe_sk_fft tsk, double parame){
  uint64_t *PK = (uint64_t *) malloc(_k * (_k + 1) * _N * sizeof(uint64_t));

  int max=_k*_k*_N;

  uniform64_distribution_vector(PK, max);

  //noise polynomial e
  uint64_t *e = (uint64_t*) malloc(_k*_N * sizeof(uint64_t));

  small_gaussian_overZ_vector(e,parame,_k*_N);

  int i, j, coeff;

  for (i = 0; i < _k*_N; ++i)
    PK[max+i] = e[i] << (64-_logBg*_ell);
  free(e);

  for (i = 0; i < _k*_k*_N; ++i)
    PK[i] &= _mask;

  double tmp;
  for (i = 0; i < _k; ++i) {
    for (j = 0; j < _k; ++j) {
      for (coeff = 0; coeff < _N; ++coeff)
        {
          fft_real[coeff] = (int64_t) PK[coeff + (j * _k + i) * _N];
        }
      memset(fft_real+_N, 0, _N * sizeof(double));
      fftw_execute(fft_forward);
      for (coeff = 0; coeff <= _N; ++coeff)
        fft_complex[coeff] *= tsk[coeff + j * (_N+1)];
      fftw_execute(fft_backwards);
      for (coeff = 0; coeff < _N; ++coeff)
        {
          tmp = fmodl(PK[max + i*_N + coeff] + (fft_real[coeff] - fft_real[_N + coeff] + _N) / (2 * _N), _two64_double);
          if (tmp < 0)
            tmp += _two64_double;
          PK[max+i*_N+coeff] = ((uint64_t) (tmp+(1L<<(63-_logBg*_ell)))) & _mask;
        }
    }
  }

  pkc_fft PK_fft = (fftw_complex*) fftw_malloc(_k * (_k + 1) * (_N + 1) * sizeof(fftw_complex));
  for (i = 0; i < _k * (_k+1); ++i) {
    for (coeff = 0; coeff < _N; ++coeff) {
      fft_real[coeff] = (int64_t) PK[coeff + i * _N];
    }
    memset(fft_real+_N, 0, _N * sizeof(double));
    fftw_execute(fft_forward);
    for (coeff = 0; coeff <= _N; ++coeff)
      PK_fft[coeff + i * (_N+1)] = fft_complex[coeff];
  }
  free(PK);
  return PK_fft;
}

void sanitize_pkc_clear_fft(pkc_fft PK) {
  fftw_free(PK);
}

void sanitize_pkc_enc_fft(tlwe_sample ct, pkc_fft PK, gaussian_param_t paramr, gaussian_param_t paramep, gaussian_param_t paramepp){
  int max = _k *_N;

  //noise polynomials
  uint64_t *e = (uint64_t*) malloc(_N * sizeof(uint64_t));
  uint64_t *r = (uint64_t*) malloc(_k * _N * sizeof(uint64_t));

  gaussian_overZ_vector(e,paramepp,_N);
  gaussian_overZ_vector(r,paramep,_k * _N);

  int i, j, coeff;

  for (i = 0; i < _k * _N; ++i)
    ct[i] += r[i] << (64-_logBg*_ell);
  for (i = 0; i < _N; ++i)
    ct[max+i] += e[i] << (64-_logBg*_ell);
  free(e);


  gaussian_overZ_vector(r,paramr,_k * _N);

  fftw_complex *r_fft = (fftw_complex*) fftw_malloc(_k * (_N + 1) * sizeof(fftw_complex));
  for (i = 0; i < _k; ++i) {
    for (coeff = 0; coeff < _N; ++coeff) {
      fft_real[coeff] = (int64_t) r[coeff + i * _N];
    }

    memset(fft_real+_N, 0, _N * sizeof(double));
    fftw_execute(fft_forward);
    for (coeff = 0; coeff <= _N; ++coeff)
      r_fft[coeff + i * (_N+1)] = fft_complex[coeff];
  }
  free(r);

  double tmp;
  for (i = 0; i < _k + 1; ++i) {
    for (coeff = 0; coeff <= _N; ++coeff)
      {
        fft_complex[coeff] = 0;
      }
    for (j = 0; j < _k; ++j) {
      for (coeff = 0; coeff <= _N; ++coeff)
        fft_complex[coeff] += PK[coeff + (i * _k + j) * (_N+1)] * r_fft[coeff + j * (_N+1)];
    }

    fftw_execute(fft_backwards);
    for (coeff = 0; coeff < _N; ++coeff)
      {
        tmp = fmodl(ct[i*_N + coeff] + (fft_real[coeff] - fft_real[_N + coeff] + _N) / (2 * _N), _two64_double);
        if (tmp < 0)
          tmp += _two64_double;
        ct[i*_N+coeff] = ((uint64_t) (tmp+(1L<<(63-_logBg*_ell)))) & _mask;
      }
  }
  free(r_fft);
}

void tgsw_encrypt_over_fft(tgsw_sample_fft ct, tlwe_sk_fft sk, double param, int* m)
{
  tgsw_sample nofft;
  tgsw_sample_init(&nofft);

  for (int i = 0; i < (_k+1) * _ell; ++i)
    tlwe_encrypt_zero_over_fft(nofft + i * ( _k + 1 ) * _N, sk, param);

  uint64_t tmp = (uint64_t) 1L << (64 - _logBg);

  for (int power = 0; power < _ell; ++power){
    for (int i = 0; i <= _k; ++i)
      for (int coeff = 0; coeff < _N; coeff++)
        nofft[power * (_k+1) * _N + i * (_ell * (_k+1) + 1) * _N + coeff] += m[coeff] * tmp;
    tmp >>= _logBg;
  }

  for (int i = 0; i < (_k+1) * (_k+1) * _ell; ++i){
    for (int coeff = 0; coeff < _N; ++coeff)
      fft_real[coeff] = nofft[i * _N + coeff] >> (64 - _split); // High order bits
    memset(fft_real+_N, 0, _N * sizeof(double));
    fftw_execute(fft_forward);
    memcpy(ct + 3 * i * (_N+1), fft_complex, (_N+1) * sizeof(fftw_complex));
    for (int coeff = 0; coeff < _N; ++coeff)
      fft_real[coeff] = (nofft[i * _N + coeff] << _split) >> (64 - _split); // Mid order bits
    memset(fft_real+_N, 0, _N * sizeof(double));
    fftw_execute(fft_forward);
    memcpy(ct + (3 * i + 1) * (_N+1), fft_complex, (_N+1) * sizeof(fftw_complex));
    for (int coeff = 0; coeff < _N; ++coeff)
      fft_real[coeff] = (nofft[i * _N + coeff] << (2 * _split)) >> (64 - _split); // Low order bits
    memset(fft_real+_N, 0, _N * sizeof(double));
    fftw_execute(fft_forward);
    memcpy(ct + (3 * i + 2) * (_N+1), fft_complex, (_N+1) * sizeof(fftw_complex));
  }

  tgsw_sample_clear(nofft);
}
tgsw_sample_fft tgsw_encrypt_fft(tlwe_sk_fft sk, double param, int* m)
{
  tgsw_sample_fft ct;
  tgsw_sample_init_fft(&ct);

  tgsw_encrypt_over_fft(ct, sk, param, m);

  return ct;
}

int* tgsw_decrypt_fft(tlwe_sk_fft sk, tgsw_sample_fft ct)
{
  tlwe_sample tlwe;
  tlwe_sample_init(&tlwe);
  for (int i = 0; i < (_k+1); ++i){
    memcpy(fft_complex, ct + 3 * (_k+1) * (_N+1) * _ell * _k + 3 * i * (_N+1), (_N+1) * sizeof(fftw_complex));
    fftw_execute(fft_backwards);
    for (int j = 0; j < _N; ++j){
      tlwe[j + i * _N] = (fft_real[j] - fft_real[_N + j] + _N) / (2 * _N);
      tlwe[j + i * _N] <<= (64 - _split);
    }
    memcpy(fft_complex, ct + 3 * (_k+1) * (_N+1) * _ell * _k + (3 * i + 1) * (_N+1), (_N+1) * sizeof(fftw_complex));
    fftw_execute(fft_backwards);
    for (int j = 0; j < _N; ++j){
      uint64_t tmp = (fft_real[j] - fft_real[_N + j] + _N) / (2 * _N);
      tlwe[j + i * _N] += tmp << (64 - 2 * _split);
    }
    memcpy(fft_complex, ct + 3 * (_k+1) * (_N+1) * _ell * _k + (3 * i + 2) * (_N+1), (_N+1) * sizeof(fftw_complex));
    fftw_execute(fft_backwards);
    for (int j = 0; j < _N; ++j){
      uint64_t tmp = (fft_real[j] - fft_real[_N + j] + _N) / (2 * _N);
      tlwe[j + i * _N] += tmp << (64 - 3 * _split);
    }
  }
  tgsw_sample_clear_fft(ct);
  int *m = tlwe_decrypt_fft(sk, _Bg, tlwe);
  return m;
}
void tgsw_decrypt_over_and_keep_fft(int *m, tlwe_sk_fft sk, tgsw_sample_fft ct)
{
  tlwe_sample tlwe;
  tlwe_sample_init(&tlwe);
  for (int i = 0; i < (_k+1); ++i){
    memcpy(fft_complex, ct + 3 * (_k+1) * (_N+1) * _ell * _k + 3 * i * (_N+1), (_N+1) * sizeof(fftw_complex));
    fftw_execute(fft_backwards);
    for (int j = 0; j < _N; ++j){
      tlwe[j + i * _N] = (fft_real[j] - fft_real[_N + j] + _N) / (2 * _N);
      tlwe[j + i * _N] <<= (64 - _split);
    }
    memcpy(fft_complex, ct + 3 * (_k+1) * (_N+1) * _ell * _k + (3 * i + 1) * (_N+1), (_N+1) * sizeof(fftw_complex));
    fftw_execute(fft_backwards);
    for (int j = 0; j < _N; ++j){
      uint64_t tmp = (fft_real[j] - fft_real[_N + j] + _N) / (2 * _N);
      tlwe[j + i * _N] += tmp << (64 - 2 * _split);
    }
    memcpy(fft_complex, ct + 3 * (_k+1) * (_N+1) * _ell * _k + (3 * i + 2) * (_N+1), (_N+1) * sizeof(fftw_complex));
    fftw_execute(fft_backwards);
    for (int j = 0; j < _N; ++j){
      uint64_t tmp = (fft_real[j] - fft_real[_N + j] + _N) / (2 * _N);
      tlwe[j + i * _N] += tmp << (64 - 3 * _split);
    }
  }
  int *res = tlwe_decrypt_fft(sk, _Bg, tlwe);
  for (int i = 0; i < _N; ++i){
    m[i] = res[i];
  }
  free(res);
}

bsk_fft generate_bsk_fft(lwe_sk sk_in, tlwe_sk_fft sk_out, double param, int n)
{
  bsk_fft BSK = (tgsw_sample_fft *) malloc(n * sizeof(tgsw_sample_fft));
  int* p = (int *) malloc(_N * sizeof(int));
  int i;
  for (i = 1; i < _N; ++i)
    p[i] = 0;

  for (i = 0; i < n; ++i)
    {
      p[0] = sk_in[i];
      BSK[i] = tgsw_encrypt_fft(sk_out, param, p);
    }
  free(p);

  return BSK;
}
void bsk_clear_fft(bsk_fft bsk, int n)
{
  for (int i = 0; i < n; ++i)
    tgsw_sample_clear_fft(bsk[i]);
  free(bsk);
}

tlwe_sample randomized_external_product_fft(tlwe_sample tlwe, tgsw_sample_fft tgsw)
{
  tlwe_sample out;
  tlwe_sample_init(&out);
  uint64_t* decomposed = (uint64_t *) malloc((_k+1) * _ell * _N * sizeof(uint64_t));

  gaussian_ginv_poly_vector(decomposed, tlwe, _k+1);

  for(int i = 0; i < (_k+1)*_N; ++i)
    out[i] = 0;
  for(int i = 0; i < (_k+1); ++i) // Column number in the matrix
    for (int j = 0; j<(_k+1)*_ell; ++j) // Row number in the matrix
      multiply_accumulate_poly_fft(out + i * _N, decomposed + j * _N, tgsw + 3 * (i + j * (_k+1)) * (_N+1));

  free(decomposed);
  return out;
}
void randomized_accum_external_product_fft(tlwe_sample out, tlwe_sample tlwe, tgsw_sample_fft tgsw)
{
  uint64_t* decomposed = (uint64_t *) malloc((_k+1) * _ell * _N * sizeof(uint64_t));

  gaussian_ginv_poly_vector(decomposed, tlwe, _k+1);

  for(int i = 0; i < (_k+1); ++i) // Column number in the matrix
    for (int j = 0; j<(_k+1)*_ell; ++j) // Row number in the matrix
      multiply_accumulate_poly_fft(out + i * _N, decomposed + j * _N, tgsw + 3 * (i + j * (_k+1)) * (_N+1));

  free(decomposed);
}
void randomized_accum_external_product_aux_fft(tlwe_sample out, tlwe_sample tlwe, tgsw_sample_fft tgsw, uint64_t* decomposed)
{
  gaussian_ginv_poly_vector(decomposed, tlwe, _k+1);

  for(int i = 0; i < (_k+1); ++i) // Column number in the matrix
    for (int j = 0; j<(_k+1)*_ell; ++j) // Row number in the matrix
      multiply_accumulate_poly_fft(out + i * _N, decomposed + j * _N, tgsw + 3 * (i + j * (_k+1)) * (_N+1));
}
void accum_external_product_fft(tlwe_sample out, tlwe_sample tlwe, tgsw_sample_fft tgsw)
{
  uint64_t* decomposed = (uint64_t *) malloc((_k+1) * _ell * _N * sizeof(uint64_t));
  decompose_tlwe(decomposed, tlwe);

  for(int i = 0; i < (_k+1); ++i) // Column number in the matrix
    for (int j = 0; j< (_k+1)*_ell; ++j) // Row number in the matrix
      multiply_accumulate_poly_fft(out + i * _N, decomposed + j * _N, tgsw + 3 * (i + j * (_k+1)) * (_N+1));

  free(decomposed);
}
tlwe_sample blind_rotate_fft(lwe_sample lwe, int n, uint64_t* testv, bsk_fft bsk)
{
  int max=_k*_N;

  int i, j, pow;

  tlwe_sample out, tmp_acc;
  tlwe_sample_init(&out);
  tlwe_sample_init(&tmp_acc);

  int bbar = (lwe[n] >> (64 - _log2N));
  for (i = 0; i < max; ++i)
    out[i] = 0;

  if (bbar < _N)
    {
      for (i = 0; i < bbar; ++i)
        out[max+i] = -testv[_N+i-bbar];
      for (; i < _N; ++i)
        out[max+i] = testv[i-bbar];
    }
  else
    {
      for (i = 0; i < bbar - _N; ++i)
        out[max+i] = testv[2*_N-bbar+i];
      for (; i < _N; ++i)
        out[max+i] = -testv[_N-bbar+i];
    }

  int abari;
  for (i = 0; i < n;  ++i)
    {
      abari = (-lwe[i] >> (64 - _log2N));
      if (abari < _N)
        {
          for (pow = 0; pow < abari; ++pow)
            for (j = 0; j <= _k; ++j)
              tmp_acc[j * _N + pow] = - out[j * _N + pow] - out[j * _N + _N + pow - abari];
          for (; pow < _N; ++pow)
            for (j = 0; j <= _k; ++j)
              tmp_acc[j * _N + pow] = - out[j * _N + pow] + out[j * _N + pow - abari];
        }
      else
        {
          for (pow = 0; pow < abari - _N; ++pow)
            for (j = 0; j <= _k; ++j)
              tmp_acc[j * _N + pow] = - out[j * _N + pow] + out[j * _N + 2*_N + pow - abari];
          for (; pow < _N; ++pow)
            for (j = 0; j <= _k; ++j)
              tmp_acc[j * _N + pow] = - out[j * _N + pow] - out[j * _N + _N + pow - abari];
        }
      accum_external_product_fft(out, tmp_acc, bsk[i]);
      
      for (j = 0; j < (_k+1)*_N; ++j)
        out[j] = (out[j]+(1L<<(63-_logBg*_ell))) & _mask;
    }


  tlwe_sample_clear(tmp_acc);

  return out;
}
tlwe_sample cp_blind_rotate_fft(lwe_sample lwe, gaussian_param_t param, int n, uint64_t* testv, bsk_fft bsk)
{
  int max=_k*_N;

  int i, j, pow;

  uint64_t decomposed[(_k + 1) * _ell * _N];
  tlwe_sample out, tmp_acc;
  tlwe_sample_init(&out);
  tlwe_sample_init(&tmp_acc);

  int bbar = (lwe[n] >> (64 - _log2N));
  for (i = 0; i < max; ++i)
    out[i] = 0;

  if (bbar < _N)
    {
      for (i = 0; i < bbar; ++i)
        out[max+i] = -testv[_N+i-bbar];
      for (; i < _N; ++i)
        out[max+i] = testv[i-bbar];
    }
  else
    {
      for (i = 0; i < bbar - _N; ++i)
        out[max+i] = testv[2*_N-bbar+i];
      for (; i < _N; ++i)
        out[max+i] = -testv[_N-bbar+i];
    }

  int abari;
  for (i = 0; i < n;  ++i)
    {
      //printf("%u\n",i);
      abari = (-lwe[i] >> (64 - _log2N));
      if (abari < _N)
        {
          for (pow = 0; pow < abari; ++pow)
            for (j = 0; j <= _k; ++j)
              tmp_acc[j * _N + pow] = - out[j * _N + pow] - out[j * _N + _N + pow - abari];
          for (; pow < _N; ++pow)
            for (j = 0; j <= _k; ++j)
              tmp_acc[j * _N + pow] = - out[j * _N + pow] + out[j * _N + pow - abari];
        }
      else
        {
          for (pow = 0; pow < abari - _N; ++pow)
            for (j = 0; j <= _k; ++j)
              tmp_acc[j * _N + pow] = - out[j * _N + pow] + out[j * _N + 2*_N + pow - abari];
          for (; pow < _N; ++pow)
            for (j = 0; j <= _k; ++j)
              tmp_acc[j * _N + pow] = - out[j * _N + pow] - out[j * _N + _N + pow - abari];
        }
      randomized_accum_external_product_aux_fft(out, tmp_acc, bsk[i], decomposed);

      for (j = 0; j < (_k+1)*_N; ++j)
        out[j] = (out[j]+(1L<<(63-_logBg*_ell))) & _mask;
    }

  uint64_t* y = (uint64_t *) malloc(_N*sizeof(uint64_t));
  gaussian_overZ_vector(y,param,_N);
  for (j = 0; j < _N; ++j)
    out[max+j] += (y[j] << (64-_logBg*_ell));
  free(y);
  tlwe_sample_clear(tmp_acc);

  return out;
}

void bootstrap_fft(lwe_sample ct, bsk_fft bsk, ksk ksk, int M_out, int n)
{
  uint64_t* test_vector=(uint64_t *) malloc(_N * sizeof(uint64_t));

  for (int i = 0; i < _N; ++i)
    test_vector[i] = _two64_double - _two64_double / (2 * M_out);


  printf("BLINDROTATE\n");

  tlwe_sample tmp_tlwe = blind_rotate_fft(ct, n, test_vector, bsk);

    //for (int i = 0; i < _N; ++i){
    //printf("%llu\n",test_vector[i]);
  //}

  tmp_tlwe[_k * _N] = fmodl(tmp_tlwe[_k * _N] + _two64_double / (2 * M_out),_two64_double); 
  //tmp_tlwe[_k * _N] += _two64_double / (2 * M_out);

  printf("EXTRACT\n");
  lwe_sample xtracted_lwe = tlwe_extract(tmp_tlwe);

  printf("KEYSWITCH\n");
  keyswitch_over_and_keep(ct, ksk,_k*_N,n,xtracted_lwe);
  tlwe_sample_clear(xtracted_lwe);

  printf("DONE\n");
  free(test_vector);
}

void sanitize_c_fft(lwe_sample ct, bsk_fft bsk, ksk ksk, pkc_fft PK, gaussian_param_t paramy, gaussian_param_t paramr, gaussian_param_t paramep, gaussian_param_t paramepp, int M_out, int n)
{
  uint64_t* test_vector=(uint64_t *) malloc(_N * sizeof(uint64_t));
  for (int i = 0; i < _N; ++i)
    test_vector[i] = _two64_double - _two64_double / (2 * M_out);

  printf("BLINDROTATE\n");
  tlwe_sample tmp_tlwe = cp_blind_rotate_fft(ct, paramy, n, test_vector, bsk);
  //
  sanitize_pkc_enc_fft(tmp_tlwe, PK, paramr, paramep, paramepp);

  tmp_tlwe[_k * _N] = fmodl(tmp_tlwe[_k * _N] + _two64_double / (2 * M_out),_two64_double); 
  //tmp_tlwe[_k * _N] += _two64_double / (2 * M_out);

  printf("EXTRACT\n");
  lwe_sample xtracted_lwe = tlwe_extract(tmp_tlwe);

  printf("KEYSWITCH\n");
  keyswitch_over_and_keep(ct, ksk,_k*_N,n,xtracted_lwe);
  tlwe_sample_clear(xtracted_lwe);

  printf("DONE\n");
  free(test_vector);
}
