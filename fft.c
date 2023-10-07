#include "fft.h"


#define _two64_double powl(2.0L,64.0L)

long double *fft_real;
fftwl_complex *fft_complex;
fftwl_plan fft_forward;
fftwl_plan fft_backwards;

void fft_init()
{
    fft_real = (long double *) fftwl_malloc(sizeof(long double) * 2 * _N);

    fft_complex = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * (_N+1));

    fftwl_import_wisdom_from_filename("wisdom");
    fft_forward = fftwl_plan_dft_r2c_1d(2 * _N, fft_real, fft_complex, FFTW_EXHAUSTIVE | FFTW_DESTROY_INPUT);
    fft_backwards = fftwl_plan_dft_c2r_1d(2 * _N, fft_complex, fft_real, FFTW_EXHAUSTIVE | FFTW_DESTROY_INPUT);
    fftwl_export_wisdom_to_filename("wisdom");
}
void fft_clear()
{
	fftwl_destroy_plan(fft_forward); fftwl_destroy_plan(fft_backwards);
    fftwl_free(fft_real); fftwl_free(fft_complex);
}

void tlwe_sk_init_fft(tlwe_sk_fft *sk)
{
	*sk = (tlwe_sk_fft) fftwl_malloc(_k * (_N+1) * sizeof(fftwl_complex));
}
tlwe_sk_fft tlwe_keygen_fft(long double param)
{
	//tlwe_sk_init(&secretkey);
	long double *key = (long double *) malloc(_N * sizeof(long double));
	tlwe_sk_fft tsk;
	tlwe_sk_init_fft(&tsk);
	size_t i, j;
	for (i = 0; i < _k; ++i)
	{
		noise_vector(key, param, _N);
		for (j = 0; j < _N; ++j)
			fft_real[j] = roundl(key[j]);
		memset(fft_real+_N, 0, _N * sizeof(long double));
		fftwl_execute(fft_forward);
		memcpy(tsk + (_N+1) * i, fft_complex, (_N+1) * sizeof(fftwl_complex));
		//for (j = 0; j < _N; ++j)
		//	secretkey[i * _N + j] = key[j];
	}
	free(key);
	return tsk;
}
void tlwe_sk_clear_fft(tlwe_sk_fft sk)
{
	fftwl_free(sk);
}

void tgsw_sample_init_fft(tgsw_sample_fft *ct)
{
	*ct = (tgsw_sample_fft) fftwl_malloc((_k+1) * (_k+1) * (_N+1) * _ell * sizeof(fftwl_complex));
}
void tgsw_sample_clear_fft(tgsw_sample_fft ct)
{
	fftwl_free(ct);
}

void multiply_accumulate_poly_fft(uint64_t *out, uint64_t *a, fftwl_complex *b)
{
	size_t i;
	for (i = 0; i < _N; ++i)
		fft_real[i] = (int64_t) a[i];
	memset(fft_real+_N, 0, _N * sizeof(long double));
	fftwl_execute(fft_forward);
	for (i = 0; i <= _N; ++i)
		fft_complex[i] *= b[i];
	fftwl_execute(fft_backwards);
	long double tmp;
	for (i = 0; i < _N; ++i)
	{
		tmp = fmodl(out[i] + (fft_real[i] - fft_real[_N + i] + _N) / (2 * _N), _two64_double);
		if (tmp < 0)
			tmp += _two64_double;
		out[i] = tmp;
	}
}

void tlwe_dotproduct_fft(tlwe_sample ct, tlwe_sk_fft sk){
	size_t i, j;
	int max=_k*_N;
	long double tmp;
	for (i = 0; i < _k; ++i)
	{
		for (j = 0; j < _N; ++j)
		{
			fft_real[j] = (int64_t) ct[j + i * _N];
		}
		memset(fft_real+_N, 0, _N * sizeof(long double));
		fftwl_execute(fft_forward);
		for (j = 0; j <= _N; ++j)
			fft_complex[j] *= sk[j + i * (_N+1)];
		fftwl_execute(fft_backwards);
		for (j = 0; j < _N; ++j)
		{
			tmp = fmodl(ct[max + j] + (fft_real[j] - fft_real[_N + j] + _N) / (2 * _N), _two64_double);
			if (tmp < 0)
				tmp += _two64_double;
			ct[max+j] = ((uint64_t) (tmp+(1<<(63-_logBg*_ell)))) & _mask;
		}
	}
}

void tlwe_encrypt_over_fft(tlwe_sample ct, tlwe_sk_fft sk, long double param, int M, int* m)
{
	int max=_k*_N;

	//noise polynomial e
    long double *e = (long double*) malloc(_N * sizeof(long double));

    noise_vector(e,param,_N);

	size_t i;

	for (i = 0; i < _N; ++i)
		ct[max+i] = ((uint64_t) (e[i] + (_two64_double* m[i])/M)) & _mask;
	free(e);

	uniform64_distribution_vector(ct, max);
	for (i = 0; i < max; ++i)
		ct[i] &= _mask;

	tlwe_dotproduct_fft(ct, sk);

}
tlwe_sample tlwe_encrypt_fft(tlwe_sk_fft sk, long double param, int M, int* m)
{
	tlwe_sample ct;
	tlwe_sample_init(&ct);
	tlwe_encrypt_over_fft(ct, sk, param, M, m);
	return ct;
}

void tlwe_encrypt_zero_over_fft(tlwe_sample ct, tlwe_sk_fft sk, long double param)
{
	int max=_k*_N;


	//noise polynomial e
    long double *e = (long double*) malloc(_N * sizeof(long double));

    noise_vector(e,param,_N);

	size_t i;

	for (i = 0; i < _N; ++i)
		ct[max+i] = ((uint64_t) e[i]) & _mask;
	free(e);

	uniform64_distribution_vector(ct, max);
	for (i = 0; i < max; ++i)
		ct[i] &= _mask;

	tlwe_dotproduct_fft(ct, sk);
}
tlwe_sample tlwe_encrypt_zero_fft(tlwe_sk_fft sk, long double param)
{
	tlwe_sample ct;
	tlwe_sample_init(&ct);
	tlwe_encrypt_zero_over_fft(ct, sk, param);
	return ct;
}

void tlwe_encrypt_zero_gaussian_over_fft(tlwe_sample ct, tlwe_sk_fft sk, gaussian_param_t param)
{
	int max=_k*_N;

	size_t i;

        gaussian_overZ_vector(ct+max,param,_N);
 	for (i = 0; i < _N; ++i)
 		ct[max+i] <<= (64-_logBg*_ell);

	uniform64_distribution_vector(ct, max);
	for (i = 0; i < max; ++i)
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
	int max=_k*_N;

	int *m = (int *) malloc(_N * sizeof(int));

	long double tmp;
	size_t i, j;
	for (i = 0; i < _k; ++i)
	{
		for (j = 0; j < _N; ++j)
			fft_real[j] = (int64_t) ct[j + i * _N];
		memset(fft_real+_N, 0, _N * sizeof(long double));
		fftwl_execute(fft_forward);
		for (j = 0; j <= _N; ++j)
			fft_complex[j] *= sk[j + i * (_N+1)];
		fftwl_execute(fft_backwards);
		for (j = 0; j < _N; ++j)
		{
			tmp = fmodl(ct[max + j] - (fft_real[j] - fft_real[_N + j] + _N) / (2 * _N), _two64_double);
			if (tmp < 0)
				tmp += _two64_double;
			ct[max + j] = tmp;
		}
	}
	for (i = 0; i < _N; ++i)
		m[i]=(int)(((long double) M * ct[max+i])/_two64_double + 0.5);
	tlwe_sample_clear(ct);
	return m;
}
void tlwe_decrypt_over_and_keep_fft(int *m, tlwe_sk_fft sk, int M, tlwe_sample ct)
{
	int max=_k*_N;

	uint64_t *b = (uint64_t *) malloc(_N * sizeof(uint64_t));
	memcpy(b, ct+max, _N * sizeof(uint64_t));

	long double tmp;
	size_t i, j;
	for (i = 0; i < _k; ++i)
	{
		for (j = 0; j < _N; ++j)
			fft_real[j] = (int64_t) ct[j + i * _N];
		memset(fft_real+_N, 0, _N * sizeof(long double));
		fftwl_execute(fft_forward);
		for (j = 0; j <= _N; ++j)
			fft_complex[j] *= sk[j + i * (_N+1)];
		fftwl_execute(fft_backwards);
		for (j = 0; j < _N; ++j)
		{
			tmp = fmodl(ct[max + j] - (fft_real[j] - fft_real[_N + j] + _N) / (2 * _N), _two64_double);
			if (tmp < 0)
				tmp += _two64_double;
			b[j] = tmp;
		}
	}
	for (i = 0; i < _N; ++i)
		m[i]=(int)(((long double) M * b[i])/_two64_double + 0.5);
	free(b);
}

lwe_sk tlwe_key_extract_fft(tlwe_sk_fft in)
{
	tlwe_sk nofft;
	tlwe_sk_init(&nofft);
	size_t i, j;
	for (i = 0; i < _k; ++i)
	{
		memcpy(fft_complex, in + (_N+1) * i, (_N+1) * sizeof(fftwl_complex));
		fftwl_execute(fft_backwards);
		for (j = 0; j < _N; ++j)
			nofft[j] = fmodl(roundl((fft_real[j] - fft_real[_N + j])/ (2 * _N)),_two64_double);
	}
	lwe_sk out = tlwe_key_extract(nofft);
	tlwe_sk_clear(nofft);
	return out;
}

/* Statistically private randomizer */
/*pkg sanitize_pk_gen_fft(tlwe_sk_fft tsk, long double param)
{
	pkg PK = (tlwe_sample *) malloc(_m * sizeof(tlwe_sample));
	size_t i;
	for (i = 0; i < _m; ++i)
    	PK[i] = tlwe_encrypt_zero_fft(tsk, param);
	return PK;
}*/

pkc_fft sanitize_pkc_gen_fft(tlwe_sk_fft tsk, long double parame){
	uint64_t *PK = (uint64_t *) malloc(_k * (_k + 1) * _N * sizeof(uint64_t));

	int max=_k*_k*_N;

	uniform64_distribution_vector(PK, max);

	//noise polynomial e
        uint64_t *e = (uint64_t*) malloc(_k*_N * sizeof(uint64_t));

        small_gaussian_overZ_vector(e,parame,_k*_N);

	size_t i, j, coeff;

        for (i = 0; i < _k*_N; ++i)
		PK[max+i] = e[i] << (64-_logBg*_ell);
	free(e);

	for (i = 0; i < _k*_k*_N; ++i)
		PK[i] &= _mask;

	long double tmp;
	for (i = 0; i < _k; ++i) {
		for (j = 0; j < _k; ++j) {
			for (coeff = 0; coeff < _N; ++coeff)
			{
				fft_real[coeff] = (int64_t) PK[coeff + (j * _k + i) * _N];
			}
			memset(fft_real+_N, 0, _N * sizeof(long double));
			fftwl_execute(fft_forward);
			for (coeff = 0; coeff <= _N; ++coeff)
				fft_complex[coeff] *= tsk[coeff + j * (_N+1)];
			fftwl_execute(fft_backwards);
			for (coeff = 0; coeff < _N; ++coeff)
			{
				tmp = fmodl(PK[max + i*_N + coeff] + (fft_real[coeff] - fft_real[_N + coeff] + _N) / (2 * _N), _two64_double);
				if (tmp < 0)
					tmp += _two64_double;
				PK[max+i*_N+coeff] = ((uint64_t) (tmp+(1<<(63-_logBg*_ell)))) & _mask;
			}
		}
	}

	pkc_fft PK_fft = (fftwl_complex*) fftwl_malloc(_k * (_k + 1) * (_N + 1) * sizeof(fftwl_complex));
	for (i = 0; i < _k * (_k+1); ++i) {
		for (coeff = 0; coeff < _N; ++coeff) {
			fft_real[coeff] = (int64_t) PK[coeff + i * _N];
		}
		memset(fft_real+_N, 0, _N * sizeof(long double));
		fftwl_execute(fft_forward);
		for (coeff = 0; coeff <= _N; ++coeff)
			PK_fft[coeff + i * (_N+1)] = fft_complex[coeff];
	}
	free(PK);
	return PK_fft;
}

void sanitize_pkc_clear_fft(pkc_fft PK) {
	fftwl_free(PK);
}

void sanitize_pkc_enc_fft(tlwe_sample ct, pkc_fft PK, gaussian_param_t paramr, gaussian_param_t paramep, gaussian_param_t paramepp){
  int max = _k *_N;

  //noise polynomials
  uint64_t *e = (uint64_t*) malloc(_N * sizeof(uint64_t));
  uint64_t *r = (uint64_t*) malloc(_k * _N * sizeof(uint64_t));

  gaussian_overZ_vector(e,paramepp,_N);
  gaussian_overZ_vector(r,paramep,_k * _N);

  size_t i, j, coeff;

  for (i = 0; i < _k * _N; ++i)
    ct[i] += r[i] << (64-_logBg*_ell);
  for (i = 0; i < _N; ++i)
    ct[max+i] += e[i] << (64-_logBg*_ell);
  free(e);


  gaussian_overZ_vector(r,paramr,_k * _N);

  fftwl_complex *r_fft = (fftwl_complex*) fftwl_malloc(_k * (_N + 1) * sizeof(fftwl_complex));
  for (i = 0; i < _k; ++i) {
    for (coeff = 0; coeff < _N; ++coeff) {
      fft_real[coeff] = (int64_t) r[coeff + i * _N];
    }

    memset(fft_real+_N, 0, _N * sizeof(long double));
    fftwl_execute(fft_forward);
    for (coeff = 0; coeff <= _N; ++coeff)
      r_fft[coeff + i * (_N+1)] = fft_complex[coeff];
  }
  free(r);

  long double tmp;
  for (i = 0; i < _k + 1; ++i) {
    for (coeff = 0; coeff <= _N; ++coeff)
      {
        fft_complex[coeff] = 0;
      }
    for (j = 0; j < _k; ++j) {
      for (coeff = 0; coeff <= _N; ++coeff)
        fft_complex[coeff] += PK[coeff + (i * _k + j) * (_N+1)] * r_fft[coeff + j * (_N+1)];
    }

    fftwl_execute(fft_backwards);
    for (coeff = 0; coeff < _N; ++coeff)
      {
        tmp = fmodl(ct[i*_N + coeff] + (fft_real[coeff] - fft_real[_N + coeff] + _N) / (2 * _N), _two64_double);
        if (tmp < 0)
          tmp += _two64_double;
        ct[i*_N+coeff] = ((uint64_t) (tmp+(1<<(63-_logBg*_ell)))) & _mask;
      }
  }
  free(r_fft);
}


void sanitize_pk_enc_online_fft(tlwe_sample out, tlwe_sk_fft tsk, long double param)
{
	size_t i, j, z;
	uint64_t r, pow;
	tlwe_sample ct;
	tlwe_sample_init(&ct);
	for (i = 0; i < _m; ++i)
	{
		random_binary(&r);
		if (r)
		{
			tlwe_encrypt_zero_over_fft(ct, tsk, param);
			uniform64_distribution(&pow);
			pow = pow >> (64 - _log2N);

			if (pow < _N)
			{
				for (z = 0; z < pow; ++z)
					for (j = 0; j <= _k; ++j)
						out[j * _N + z] -= ct[j * _N + _N + z - pow];
				for (; z < _N; ++z)
					for (j = 0; j <= _k; ++j)
						out[j * _N + z] += ct[j * _N + z - pow];
			}
			else
			{
				for (z = 0; z < pow - _N; ++z)
					for (j = 0; j <= _k; ++j)
						out[j * _N + z] += ct[j * _N + 2*_N + z - pow];
				for (; z < _N; ++z)
					for (j = 0; j <= _k; ++j)
						out[j * _N + z] -= ct[j * _N + _N + z - pow];
			}
		}
	}
    free(ct);
}


void tgsw_encrypt_over_fft(tgsw_sample_fft ct, tlwe_sk_fft sk, long double param, int* m)
{
	tgsw_sample nofft;
	tgsw_sample_init(&nofft);

	size_t i;
	for (i = 0; i < (_k+1) * _ell; ++i)
		tlwe_encrypt_zero_over_fft(nofft + i * ( _k + 1 ) * _N, sk, param);

	uint64_t tmp = (uint64_t) 1 << (64 - _logBg);

	size_t coeff;
	for (size_t power = 0; power < _ell; ++power)
	{
		for (i = 0; i <= _k; ++i)
		{
			for (coeff = 0; coeff < _N; coeff++)
			nofft[power * (_k+1) * _N + i * (_ell * (_k+1) + 1) * _N + coeff] += m[coeff] * tmp;
		}
		tmp >>= _logBg;
	}

	for (i = 0; i < (_k+1) * (_k+1) * _ell; ++i)
	{
		for (coeff = 0; coeff < _N; ++coeff)
			fft_real[coeff] = (int64_t) nofft[i * _N + coeff];
		memset(fft_real+_N, 0, _N * sizeof(long double));
		fftwl_execute(fft_forward);
		memcpy(ct + i * (_N+1), fft_complex, (_N+1) * sizeof(fftwl_complex));
	}

	tgsw_sample_clear(nofft);
}
tgsw_sample_fft tgsw_encrypt_fft(tlwe_sk_fft sk, long double param, int* m)
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
	size_t i, j;
	for (i = 0; i < (_k+1); ++i)
	{
		memcpy(fft_complex, ct + (_k+1) * (_N+1) * _ell * _k + i * (_N+1), (_N+1) * sizeof(fftwl_complex));
		fftwl_execute(fft_backwards);
		for (j = 0; j < _N; ++j)
			tlwe[j + i * _N] = fmodl((fft_real[j] - fft_real[_N + j] + _N) / (2 * _N),_two64_double);
	}
	tgsw_sample_clear_fft(ct);
	int *m = tlwe_decrypt_fft(sk, _Bg, tlwe);
	return m;
}
void tgsw_decrypt_over_and_keep_fft(int *m, tlwe_sk_fft sk, tgsw_sample_fft ct)
{
	tlwe_sample tlwe;
	tlwe_sample_init(&tlwe);
	size_t i, j;
	for (i = 0; i < (_k+1); ++i)
	{
		memcpy(fft_complex, ct + (_k+1) * (_N+1) * _ell * _k + i * (_N+1), (_N+1) * sizeof(fftwl_complex));
		fftwl_execute(fft_backwards);
		for (j = 0; j < _N; ++j)
			tlwe[j + i * _N] = fmodl((fft_real[j] - fft_real[_N + j] + _N) / (2 * _N),_two64_double);
	}
	tlwe_decrypt_over_and_keep_fft(m, sk, _Bg, tlwe);
}

bsk_fft generate_bsk_fft(lwe_sk sk_in, tlwe_sk_fft sk_out, long double param, size_t n)
{
	bsk_fft BSK = (tgsw_sample_fft *) malloc(n * sizeof(tgsw_sample_fft));
	int* p = (int *) malloc(_N * sizeof(int));
	size_t i;
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
void bsk_clear_fft(bsk_fft bsk, size_t n)
{
	for (size_t i = 0; i < n; ++i)
		tgsw_sample_clear_fft(bsk[i]);
	free(bsk);
}

tlwe_sample randomized_external_product_fft(tlwe_sample tlwe, tgsw_sample_fft tgsw)
{
	tlwe_sample out;
	tlwe_sample_init(&out);
	uint64_t* decomposed = (uint64_t *) malloc((_k+1) * _ell * _N * sizeof(uint64_t));
	size_t i, j;

	gaussian_ginv_poly_vector(decomposed, tlwe, _k+1);

	for(i = 0; i < (_k+1)*_N; ++i)
			out[i] = 0;
	for(i = 0; i < (_k+1); ++i) // Column number in the matrix
		for (j = 0; j<(_k+1)*_ell; ++j) // Row number in the matrix
			multiply_accumulate_poly_fft(out + i * _N, decomposed + j * _N, tgsw + (i + j * (_k+1)) * (_N+1));

	free(decomposed);
	return out;
}
void randomized_accum_external_product_fft(tlwe_sample out, tlwe_sample tlwe, tgsw_sample_fft tgsw)
{
	uint64_t* decomposed = (uint64_t *) malloc((_k+1) * _ell * _N * sizeof(uint64_t));
	size_t i, j;

	gaussian_ginv_poly_vector(decomposed, tlwe, _k+1);

	for(i = 0; i < (_k+1); ++i) // Column number in the matrix
		for (j = 0; j<(_k+1)*_ell; ++j) // Row number in the matrix
			multiply_accumulate_poly_fft(out + i * _N, decomposed + j * _N, tgsw + (i + j * (_k+1)) * (_N+1));

	free(decomposed);
}
void randomized_accum_external_product_aux_fft(tlwe_sample out, tlwe_sample tlwe, tgsw_sample_fft tgsw, uint64_t* decomposed)
{
	size_t i, j;

	gaussian_ginv_poly_vector(decomposed, tlwe, _k+1);

	for(i = 0; i < (_k+1); ++i) // Column number in the matrix
		for (j = 0; j<(_k+1)*_ell; ++j) // Row number in the matrix
			multiply_accumulate_poly_fft(out + i * _N, decomposed + j * _N, tgsw + (i + j * (_k+1)) * (_N+1));
}
void accum_external_product_fft(tlwe_sample out, tlwe_sample tlwe, tgsw_sample_fft tgsw)
{
	uint64_t* decomposed = (uint64_t *) malloc((_k+1) * _ell * _N * sizeof(uint64_t));
	size_t i, j;

	decompose_tlwe(decomposed, tlwe);

	for(i = 0; i < (_k+1); ++i) // Column number in the matrix
		for (j = 0; j< (_k+1)*_ell; ++j) // Row number in the matrix
			multiply_accumulate_poly_fft(out + i * _N, decomposed + j * _N, tgsw + (i + j * (_k+1)) * (_N+1));

	free(decomposed);
}
tlwe_sample blind_rotate_fft(lwe_sample lwe, int M, size_t n, uint64_t* testv, bsk_fft bsk)
{
	size_t max=_k*_N;

	size_t i, j, pow;

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
			out[j] = (out[j]+(1<<(63-_logBg*_ell))) & _mask;
	}


	//lwe_sample_clear(lwe);
	tlwe_sample_clear(tmp_acc);

	return out;
}
tlwe_sample cp_blind_rotate_fft(lwe_sample lwe, gaussian_param_t param, int M, size_t n, uint64_t* testv, bsk_fft bsk)
{
	size_t max=_k*_N;

	size_t i, j, pow;

	uint64_t* decomposed = (uint64_t *) malloc((_k+1) * _ell * _N * sizeof(uint64_t));
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
			out[j] = (out[j]+(1<<(63-_logBg*_ell))) & _mask;
	}

	uint64_t* y = (uint64_t *) malloc(_N*sizeof(uint64_t));
        gaussian_overZ_vector(y,param,_N);
        for (j = 0; j < _N; ++j)
          out[max+j] += (y[j] << (64-_logBg*_ell));
	free(y);
	tlwe_sample_clear(tmp_acc);
	free(decomposed);

	return out;
}

/*tlwe_sample cps_blind_rotate_online_fft(lwe_sample lwe, long double param, long double noiseparam, int M, size_t n, uint64_t* testv, bsk_fft bsk, tlwe_sk_fft sk)
{
	size_t max=_k*_N;

	size_t i, j, pow;

	tlwe_sample out, tmp_acc, pkenc;
	uint64_t* y = (uint64_t *) malloc(_N*sizeof(uint64_t));
	tlwe_sample_init(&out);
	tlwe_sample_init(&tmp_acc);
	tlwe_sample_init(&pkenc);

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
		randomized_accum_external_product_fft(out, tmp_acc, bsk[i], param);

		gaussian_overZ_vector(y,param,_N);
		for (j = 0; j < _N; ++j)
			out[max+j] += (y[j] << (64-_logBg*_ell));

		tlwe_encrypt_zero_over_fft(pkenc, sk, noiseparam);
		for (j = 0; j < (_k+1)*_N; ++j)
			out[j] += pkenc[j];

		for (j = 0; j < (_k+1)*_N; ++j)
			out[j] = (out[j]+(1<<(63-_logBg*_ell))) & _mask;
	}

	free(y);
	tlwe_sample_clear(tmp_acc);
	tlwe_sample_clear(pkenc);


	return out;
}*/

void bootstrap_fft(lwe_sample ct, bsk_fft bsk, ksk ksk, int M_out, size_t n)
{
	uint64_t* test_vector=(uint64_t *) malloc(_N * sizeof(uint64_t));
	size_t i;
	for (i = 0; i < _N; ++i)
		test_vector[i] = - _two64_double / (2 * M_out);

	printf("BLINDROTATE\n");
	tlwe_sample tmp_tlwe = blind_rotate_fft(ct, M_out,n,test_vector,bsk);

	tmp_tlwe[_k * _N] += _two64_double / (2 * M_out);

	printf("EXTRACT\n");
	lwe_sample xtracted_lwe = tlwe_extract(tmp_tlwe);

	printf("KEYSWITCH\n");
	keyswitch_over_and_keep(ct, ksk,_k*_N,n,xtracted_lwe);
	tlwe_sample_clear(xtracted_lwe);

	printf("DONE\n");
	free(test_vector);
}

/*void sanitize_s_fft(lwe_sample ct, bsk_fft bsk, ksk ksk, pks PK, long double param, int M_out, size_t n)
{
	uint64_t* test_vector=(uint64_t *) malloc(_N * sizeof(uint64_t));
	size_t i;
	for (i = 0; i < _N; ++i)
		test_vector[i] = - _two64_double / (2 * M_out);

	printf("BLINDROTATE\n");
	tlwe_sample tmp_tlwe = cp_blind_rotate_fft(ct, param, M_out,n,test_vector,bsk);

	tmp_tlwe[_k * _N] += _two64_double / 8;

	printf("EXTRACT\n");
	lwe_sample xtracted_lwe = tlwe_extract(tmp_tlwe);

	printf("KEYSWITCH\n");
	keyswitch_over_and_keep(ct, ksk,_k*_N,n,xtracted_lwe);
	tlwe_sample_clear(tmp_tlwe);

	printf("DONE\n");
	free(test_vector);
}*/

void sanitize_c_fft(lwe_sample ct, bsk_fft bsk, ksk ksk, pkc_fft PK, gaussian_param_t paramy, gaussian_param_t paramr, gaussian_param_t paramep, gaussian_param_t paramepp, int M_out, size_t n)
{
	uint64_t* test_vector=(uint64_t *) malloc(_N * sizeof(uint64_t));
	size_t i;
	for (i = 0; i < _N; ++i)
		test_vector[i] = - _two64_double / (2 * M_out);

	printf("BLINDROTATE\n");
	tlwe_sample tmp_tlwe = cp_blind_rotate_fft(ct, paramy, M_out,n,test_vector,bsk);
	//
	sanitize_pkc_enc_fft(tmp_tlwe, PK, paramr, paramep, paramepp);

	tmp_tlwe[_k * _N] += _two64_double / (2 * M_out);

	printf("EXTRACT\n");
	lwe_sample xtracted_lwe = tlwe_extract(tmp_tlwe);
	//tlwe_sample_clear(tmp_tlwe);

	printf("KEYSWITCH\n");
	keyswitch_over_and_keep(ct, ksk,_k*_N,n,xtracted_lwe);
	tlwe_sample_clear(xtracted_lwe);
	
	printf("DONE\n");
	free(test_vector);

	
}
