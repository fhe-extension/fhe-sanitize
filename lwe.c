/* Through this file, n denotes the dimension of the underlying LWE problem */
/* The keys are random binary vectors of size n, we note them s */
/* The ciphertexts are composed of n + 1 values */
/* The first n values of a ciphertext are uniformly random, we call that part a*/
/* The last value of a ciphertext is sa + e + m * 2^64/M */
/* Where M is a bound on the message space */
/* Memory is allocated during key generation and decryption */
/* Memory is deallocated during decryption */
/* Key switching deallocates input and allocate memory for output ciphertext */
/* Variants of KS, Enc, and Dec are proposed that do not allocate or free memory */

#include "lwe.h"

/* Value of 2^64 as a long double, used during encryption */
#define _two64_double powl(2.0L,64.0L) 


/* Functions that manage the memory usage of ciphertexts */
void lwe_sample_init(lwe_sample *ct, size_t n)
{
	*ct = (lwe_sample) malloc((n+1) * sizeof(uint64_t)); // 
}
void lwe_sample_clear(lwe_sample ct)
{
	free(ct);
}


/* Functions that manage the memory usage of secret keys */
void lwe_sk_init(lwe_sk *s, size_t n)
{
	*s = (lwe_sk) malloc (n * sizeof(uint64_t)); // Keys are size n
}
void lwe_sk_clear(lwe_sk sk)
{
	free(sk);
}


/* Function that zeroes out a ciphertext */
void lwe_sample_zero(lwe_sample ct, size_t n)
{
	for(size_t i=0; i <= n; ++i) {
		ct[i] = 0;
	}
}


/* Generates a secret key and outputs it */
lwe_sk lwe_keygen(size_t n)
{
	lwe_sk s;
    lwe_sk_init(&s, n);
	random_binary_vector(s, n);
	return s;
}


/* Generates a ciphertext for message m with message space M and error parameter param and outputs it */
lwe_sample lwe_encrypt(lwe_sk s, int M, long double param, int mu, size_t n)
{
	lwe_sample ct;
	lwe_sample_init(&ct, n);
	uniform64_distribution_vector(ct, n);
   	long double e; 
	noise(&e, param);
	ct[n] = (uint64_t) (e + (mu * _two64_double)/M);
	size_t i;
	for (i = 0; i < n; ++i)
	{
		ct[n] += ct[i]*s[i];
	}
	return ct;
}

/* Encrypt without allocating memory */
void lwe_encrypt_over(lwe_sample ct, lwe_sk s, int M, long double param, int mu, size_t n)
{	
	uniform64_distribution_vector(ct, n);
   	long double e; 
	noise(&e, param);
	ct[n] = (uint64_t) (e + (mu * _two64_double)/M);
	size_t i;
	for (i = 0; i < n; ++i)
	{
		ct[n] += ct[i]*s[i];
	}
}


/* Decrypts ciphertext ct with secret key n and message space M and outputs the result, cleans the ciphertext */
int lwe_decrypt(lwe_sk s, int M, lwe_sample ct, size_t n)
{
	uint64_t b = ct[n];
	
	size_t i;
	for (i = 0; i < n; ++i){
		b -= ct[i] * s[i];
	}

	lwe_sample_clear(ct);
    
    return (int) (((long double) b * M) / _two64_double + 0.5); 
}

/* Decrypt without freeing the ciphertext */
int lwe_decrypt_and_keep(lwe_sk s, int M, lwe_sample ct, size_t n)
{
	uint64_t b = ct[n];
	size_t i;
	for (i = 0; i < n; ++i){
		b -= ct[i] * s[i];
	}

    return (int) (((long double) b * M) / _two64_double + 0.5); 
}

/* Generates a key switching key and outputs it */
ksk generate_ksk(lwe_sk sk_in, lwe_sk sk_out, long double param, size_t n_in, size_t n_out)
{
	ksk ksk = (lwe_sample*) malloc(n_in * _t * sizeof(lwe_sample));;
	size_t i, j;
	for(i = 0; i < n_in; ++i)
		for(j = 0; j < _t; ++j)
			/* KSK_i,j = LWE(s_i * 2^64/Bks^(j+1)) */
			ksk[i*_t + j] = lwe_encrypt(sk_out, pow(_Bks, (j+1)), param, sk_in[i], n_out);
	return ksk;
}

/* Cleans up key switching key */
void ksk_clear(ksk ksk, size_t n_in, size_t n_out)
{
	size_t i, j;
	for (i = 0; i < n_in; ++i)
		for (j = 0; j < _t; ++j)
			lwe_sample_clear(ksk[i*_t + j]);
	free(ksk);
}



lwe_sample keyswitch(ksk ksk, size_t n_in, size_t n_out, lwe_sample in)
{
	const uint64_t mask = (1 << _logBks) - 1;

	lwe_sample out;
	lwe_sample_init(&out, n_out);

	size_t i,j,k;
	for(i = 0; i < n_out; ++i) {
		out[i] = 0;
	}
	out[n_out] = in[n_in];

	uint64_t aij;

	for(i = 0; i < n_in; ++i)
		for(j = 0; j < _t; ++j) {
		//decomposition of ai in base base, ..,base^{t}
			aij= (in[i]>>(64-(j+1)*_logBks)) & mask ; // aij = ai*Bks^(j+1)/q
			if(aij != 0)
				for(k = 0; k < n_out + 1; ++k)
					out[k] -= aij * (ksk[i*_t + j])[k];
		}
	lwe_sample_clear(in);
	return out;
}

void keyswitch_over_and_keep(lwe_sample out, ksk ksk, size_t n_in, size_t n_out, lwe_sample in)
{
	const uint64_t mask = (1 << _logBks) - 1;

	size_t i,j,k;
	for(i = 0; i < n_out; ++i) {
		out[i] = 0;
	}
	out[n_out] = in[n_in];

	uint64_t aij;

	for(i = 0; i < n_in; ++i)
		for(j = 0; j < _t; ++j) {
		//decomposition of ai in base base, ..,base^{t}
			aij= (in[i]>>(64-(j+1)*_logBks)) & mask ; // aij = ai*Bks^(j+1)/q
			if(aij != 0)
				for(k = 0; k < n_out + 1; ++k)
					out[k] -= aij * (ksk[i*_t + j])[k];
		}
}







