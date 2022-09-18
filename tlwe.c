//#pragma once

#include "tlwe.h"

#define _two64_double powl(2.0L,64.0L)

#define _pool_size 1//21420

//Precomputed pools of encryptions of zero
struct pkenc_pool
{
	size_t size;
	tlwe_sample* pool;
};
struct pkenc_pool Pool; // pointer pools is initialized to NULL
// Hardcoded maxSize in the init function

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
	//*ct = (tlwe_sample) malloc((k+1) * N * sizeof(uint64_t));
	size_t max = (_k+1)*_N;
	for(size_t i=0; i < max; ++i) {
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
	for (size_t i = 0; i < (_k + 1) * _N; ++i)
		out[i] += a[i];
}

// out += a * b // Inputs are degree N poly modulo X^N+1
void multiply_accumulate_poly(uint64_t *out, uint64_t *a, uint64_t *b)
{
	size_t coeff1, coeff2;
	for (coeff1 = 0; coeff1 < _N; ++coeff1)
	{
		for (coeff2 = 0; coeff2 <= coeff1; coeff2++)
	    	out[coeff1] += a[coeff2]*b[coeff1-coeff2];
		for (; coeff2 < _N; coeff2++)
	    	out[coeff1] -= a[coeff2]*b[_N+coeff1-coeff2];
	}
}

/*void tlwe_add_mul(tlwe_sample out, tlwe_sample a, uint64_t* b, size_t k, size_t N){
	size_t j;
	uint64_t *tmpa = (uint64_t *) malloc(N * sizeof(uint64_t));		

		//for( pow= 0; pow < N; ++pow) 
		//	tmpout[pow]=0; 

	for( j= 0; j <= k; ++j) {
		memcpy(tmpa, a+j, N * sizeof(uint64_t));
		multiply_accumulate_poly(out+j*N, tmpa, b, N);
	}
free(tmpa);
}*/


tlwe_sample tlwe_encrypt(tlwe_sk sk, long double param, int M, int *m)
{   
	int max=_k*_N;

	tlwe_sample ct;
	tlwe_sample_init(&ct);

	//noise polynomial e
    long double *e = (long double*) malloc(_N * sizeof(long double));

    noise_vector(e,param,_N);    

	size_t poly, coeff1, coeff2;
	
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
	return ct;
}
void tlwe_encrypt_over(tlwe_sample ct, tlwe_sk sk, long double param, int M, int *m)
{   
	//tlwe_sample_zero(ct, k, N);

	int max=_k*_N;

	 //noise polynomial e
    long double *e = (long double*) malloc(_N * sizeof(long double));

    //b de taille N b=e gaussien
    noise_vector(e,param,_N);    

	size_t poly, coeff1, coeff2;
	
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
tlwe_sample tlwe_encrypt_zero(tlwe_sk sk, long double param)
{   
	int max=_k*_N;

	tlwe_sample ct;
	tlwe_sample_init(&ct);

	//noise polynomial e
    long double *e = (long double*) malloc(_N * sizeof(long double));

    noise_vector(e,param,_N);    

	size_t poly, coeff1, coeff2;
	
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
	return ct;
}
void tlwe_encrypt_zero_over(tlwe_sample ct, tlwe_sk sk, long double param)
{   
	int max=_k*_N;

	//noise polynomial e
    long double *e = (long double*) malloc(_N * sizeof(long double));

    noise_vector(e,param,_N);    

	size_t poly, coeff1, coeff2;
	
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
tlwe_sample tlwe_encrypt_zero_gaussian(tlwe_sk sk, long double param)
{   
	int max=_k*_N;

	tlwe_sample ct;
	tlwe_sample_init(&ct);

	size_t poly, coeff1, coeff2;

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
	return ct;
}
void tlwe_encrypt_zero_gaussian_over(tlwe_sample ct, tlwe_sk sk, long double param)
{   
	int max=_k*_N;

	size_t poly, coeff1, coeff2;		

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

int* tlwe_decrypt(tlwe_sk sk, int M, tlwe_sample ct)
{
	int max=_k*_N;

	int *m = (int *) malloc(_N * sizeof(int));
	
	size_t poly, coeff1, coeff2;
	for (poly = 0; poly < _k; ++poly)
		for (coeff1 = 0; coeff1 < _N; ++coeff1)
		{
			for (coeff2 = 0; coeff2 <= coeff1; coeff2++)
		    	ct[max+coeff1] -= ct[poly*_N+coeff2]*sk[poly*_N+coeff1-coeff2];
			for (; coeff2 < _N; coeff2++)	   
		    	ct[max+coeff1] += ct[poly*_N+coeff2]*sk[poly*_N+_N+coeff1-coeff2];
		}
	for (coeff1 = 0; coeff1 < _N; ++coeff1)
		m[coeff1]=(int)(((long double) M * ct[max+coeff1])/_two64_double + 0.5);
	tlwe_sample_clear(ct);
	return m;
}

void tlwe_decrypt_over_and_keep(int* m, tlwe_sk sk, int M, tlwe_sample ct)
{
	int max=_k*_N;

	uint64_t *b = (uint64_t *) malloc(_N * sizeof(uint64_t));
	memcpy(b, ct+max, _N * sizeof(uint64_t));
	
	size_t poly, coeff1, coeff2;
	for (poly = 0; poly < _k; ++poly)
		for (coeff1 = 0; coeff1 < _N; ++coeff1)
		{
			for (coeff2 = 0; coeff2 <= coeff1; coeff2++)
		    	b[coeff1] -= ct[poly*_N+coeff2]*sk[poly*_N+coeff1-coeff2];
			for (; coeff2 < _N; coeff2++)
		    	b[coeff1] += ct[poly*_N+coeff2]*sk[poly*_N+_N+coeff1-coeff2];
		}
	for (coeff1 = 0; coeff1 < _N; ++coeff1)
		m[coeff1]=(int)(((long double) M * b[coeff1])/_two64_double + 0.5);
	free(b);
}


lwe_sample tlwe_extract(tlwe_sample in)
{
	lwe_sample out;
	lwe_sample_init(&out, _k * _N);
	out[_k*_N] = in[_k*_N];  

	size_t i,j;
	size_t a_pos = 0;
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



/*
void tlwe_sample_extract(lwe_sample *out, tlwe_sample in, size_t k, size_t N)
{
	//size=k*N to (k+1)*N
	//out->b = in.b.coeffs[0];
	for(int j = 0; j < N; ++j) 
	{
		(*out)[k*N+j-1] = (uint64_t) in[k*N+j-1]; 
	}
	
	int a_pos = 0;
	for(int i = 0; i < k; ++i) {
		//out[a_pos++] = in.a[i].coeffs[0];  // add the constant term as it is
		// reverse and change the sign of the others
		for(int j = N - 1; j > 0; --j) {
			(*out)[a_pos++] = -in[i*N+j];
		}
	}
}
*/
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
			sk_out[sk_pos++] = (uint64_t) in[i*_N+j]; 
		}
	}
	return sk_out;
}


//B and ell hard coded
//return a unint64* table of size (N*ell)*(k+1) = e_11, ..., ek+1,ell in R_q
//s.t. 
void decompose_poly(uint64_t *out, uint64_t *in)
{
	size_t coeff, i;
	for (coeff = 0; coeff < _N; ++coeff)
		for (i=0; i < _ell; ++i)
			out[i*_N+coeff] = (in[coeff] << i * _logBg) >> (64 - _logBg);
}
void decompose_tlwe(uint64_t *out, tlwe_sample in)
{
	for (size_t i = 0; i < _k+1; ++i)
		decompose_poly(out+(i*_ell*_N), in + i*_N);
}

//sanitize keys  

/*void sanitize_pkg_init(tlwe_sample* PK, size_t m, size_t k, size_t N ){
	
	PK = (tlwe_sample*) malloc (m *  sizeof(tlwe_sample));
	for (size_t i=0;i<m; i++) 
       tlwe_sample_init(PK+i,k,N);
}*/

pkg sanitize_pk_gen(tlwe_sk tsk, long double param) 
{
	pkg PK = (tlwe_sample *) malloc(_m * sizeof(tlwe_sample));
	size_t i;
	for (i = 0; i < _m; ++i)
    	PK[i] = tlwe_encrypt_zero(tsk, param);
	return PK;
}


void sanitize_pk_clear(pkg PK)
{
	for (size_t i = 0; i < _m; ++i) 
     	tlwe_sample_clear(PK[i]);
     free(PK);
}

/*void sanitize_epk_init(tlwe_sample *UPK, size_t k, size_t N ){
	*UPK = (tlwe_sample) malloc (sizeof(tlwe_sample));
}



void sanitize_epk_clear(epk UPK){
	//for (i=0;i<m; i++) 
	tlwe_sample_clear(UPK);
}*/


//Computes an encryption of zero using PK and adds it to ct
void sanitize_pk_enc_empty_pool(tlwe_sample ct, pkg PK)
{
	size_t i, j, z;
	uint64_t* r = (uint64_t*) malloc(_m * sizeof(uint64_t));
	uint64_t* pow = (uint64_t*) malloc(_m * sizeof(uint64_t));
	random_binary_vector(r, _m);
	uniform64_distribution_vector(pow, _m);
	for (i = 0; i < _m; ++i)
	{
		if (r[i])
		{
			pow[i] = pow[i] >> (64 - _log2N);
			if (pow[i] < _N)
			{
				for (z = 0; z < pow[i]; ++z)
					for (j = 0; j <= _k; ++j)
						ct[j * _N + z] -= PK[i][j * _N + _N + z - pow[i]];
				for (; z < _N; ++z)
					for (j = 0; j <= _k; ++j)
						ct[j * _N + z] += PK[i][j * _N + z - pow[i]];
			}
			else
			{
				for (z = 0; z < pow[i] - _N; ++z)
					for (j = 0; j <= _k; ++j)
						ct[j * _N + z] += PK[i][j * _N + 2*_N + z - pow[i]];
				for (; z < _N; ++z)
					for (j = 0; j <= _k; ++j)
						ct[j * _N + z] -= PK[i][j * _N + _N + z - pow[i]];
			}
		}
	}
    free(r);
    free(pow);
}
//Adds an encryption of zero to the ciphertext ct, precomputed if available
void sanitize_pk_enc(tlwe_sample ct, pkg PK)
{
	if (Pool.pool != NULL && Pool.size > 0)
	{
		//tlwe_add_to(ct, Pool.pool[--Pool.size]); // This line for real precomputation
		//tlwe_sample_clear(Pool.pool[Pool.size]); // This line for real precomputation
		tlwe_add_to(ct, Pool.pool[Pool.size-1]); // This line for dummy precomputation
		return;
	}
	sanitize_pk_enc_empty_pool(ct, PK);
}

//Precomputes _pool_size (defined at the top of this file) encryptions of zero
void precompute_pkenc(pkg PK)
{
	Pool.size = _pool_size;
	if (Pool.pool == NULL)
		Pool.pool = (tlwe_sample *) malloc(_pool_size * sizeof(tlwe_sample));
	for (size_t i = 0; i < _pool_size; ++i)
	{
		tlwe_sample_init(Pool.pool + i);
		tlwe_sample_zero(Pool.pool[i]);
		sanitize_pk_enc_empty_pool(Pool.pool[i], PK);
	}
}

void clear_pkenc()
{
	if (Pool.pool != NULL)
	{
		free(Pool.pool);
		Pool.pool = NULL;
	}
}
