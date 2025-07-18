
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

/* Value of 2^64 as a double, used during encryption */
#define _two64_double pow(2.0L,64.0L) 


/* Functions that manage the memory usage of ciphertexts */
void lwe_sample_init(lwe_sample *ct, int n)
{
  *ct = (lwe_sample) malloc((n+1) * sizeof(uint64_t)); // 
}
void lwe_sample_clear(lwe_sample ct)
{
  free(ct);
}


/* Functions that manage the memory usage of secret keys */
void lwe_sk_init(lwe_sk *s, int n)
{
  *s = (lwe_sk) malloc (n * sizeof(uint64_t)); // Keys are size n
}
void lwe_sk_clear(lwe_sk sk)
{
  free(sk);
}


/* Function that zeroes out a ciphertext */
void lwe_sample_zero(lwe_sample ct, int n)
{
  for(int i=0; i <= n; ++i) {
    ct[i] = 0;
  }
}


/* Generates a secret key and outputs it */
lwe_sk lwe_keygen(int n)
{
  lwe_sk s;
  lwe_sk_init(&s, n);
  uint64_t coins[n];
  random_binary_vector(coins, n);
  for (int i = 0; i < n; ++i)
    s[i] = coins[i];
  return s;
}

lwe_sk lwe_gaussian_keygen(int n, double param)
{
  lwe_sk s;
  lwe_sk_init(&s, n);
  double *g = malloc(n*sizeof(double));
  noise_vector(g, param, n);
  for (int i = 0; i < n; ++i)
    s[i] = g[i] < 0 ? round(g[i] + _two64_double) : round(g[i]);
  free(g);
  return s;
}

/* Encrypt without allocating memory */
void lwe_encrypt_over(lwe_sample ct, lwe_sk s, int M, double param, int mu, int n)
{
  mu = mu % M;
  if (mu < 0) mu += M;
  uniform64_distribution_vector(ct, n);
  double e;
  noise(&e, param);
  e += (mu * _two64_double)/M;
  ct[n] = e < 0 ? round(e + _two64_double) : round(e);
  for (int i = 0; i < n; ++i) {
    ct[n] += ct[i]*s[i];
  }
}

/* Generates a ciphertext for message m with message space M and error parameter param and outputs it */
lwe_sample lwe_encrypt(lwe_sk s, int M, double param, int mu, int n)
{
  lwe_sample ct;
  lwe_sample_init(&ct, n);

  lwe_encrypt_over(ct, s, M, param, mu, n);

  return ct;
}

/* Decrypts ciphertext ct with secret key n and message space M and outputs the result, cleans the ciphertext */
int lwe_decrypt(lwe_sk s, int M, lwe_sample ct, int n)
{
  uint64_t b = ct[n];

  for (int i = 0; i < n; ++i){
    b -= ct[i] * s[i];
  }

  lwe_sample_clear(ct);

  return (int) (((double) b * M) / _two64_double + 0.5);
}

/* Decrypt without freeing the ciphertext */
int lwe_decrypt_and_keep(lwe_sk s, int M, lwe_sample ct, int n, double *err)
{
  uint64_t b = ct[n];

  for (int i = 0; i < n; ++i){
    b -= ct[i] * s[i];
  }
  double m = round(((double) b * M) / _two64_double);
  *err = (((double) b) - _two64_double * (double) m / M) / _two64_double;

  return (int) (((double) b * M) / _two64_double + 0.5);
}

/* Generates a key switching key and outputs it */
ksk generate_ksk(lwe_sk sk_in, lwe_sk sk_out, double param, int n_in, int n_out)
{
  ksk ksk = (lwe_sample*) malloc(n_in * _t * sizeof(lwe_sample));;
  for(int i = 0; i < n_in; ++i)
    for(int j = 0; j < _t; ++j)
      /* KSK_i,j = LWE(s_i * 2^64/Bks^(j+1)) */
      ksk[i*_t + _t - j - 1] = lwe_encrypt(sk_out, pow(_Bks, (j+1)), param, sk_in[i], n_out);
  return ksk;
}

/* Cleans up key switching key */
void ksk_clear(ksk ksk, int n_in)
{
  for (int i = 0; i < n_in; ++i)
    for (int j = 0; j < _t; ++j)
      lwe_sample_clear(ksk[i*_t + j]);
  free(ksk);
}



lwe_sample keyswitch(ksk ksk, int n_in, int n_out, lwe_sample in)
{
  lwe_sample out;
  lwe_sample_init(&out, n_out);

  keyswitch_over_and_keep(out, ksk, n_in, n_out, in);

  lwe_sample_clear(in);
  return out;
}

void keyswitch_over_and_keep(lwe_sample out, ksk ksk, int n_in, int n_out, lwe_sample in)
{
  //const uint64_t mask = (1 << _logBks) - 1;

  for(int i = 0; i < n_out; ++i) {
    out[i] = 0;
  }
  out[n_out] = in[n_in];

  uint64_t ai;
  uint64_t aij;
  //int64_t approx;

  for(int i = 0; i < n_in; ++i)  {
    ai = in[i] + (1L << (63 - _t * _logBks));
    ai >>= (64 - _t * _logBks);
    //printf("ai = %lld\n", in[i]);
    for(int j = 0; j < _t; ++j) {
      //decomposition of ai in base base, ..,base^{t}
      int64_t tmp = ai << (64 - _logBks);
      aij = tmp >> (64 - _logBks); // j=0 LSB j=_t-1 MSB
      //printf("aij = %d, i = %d, j=%d, \n", aij, i, j);
      if(aij != 0) {
        //approx += aij * (1L << (64 - _t * _logBks + j));
        for(int k = 0; k < n_out + 1; ++k)
          out[k] -= aij * (ksk[i*_t + j])[k];
      }
      ai -= aij;
      ai >>= _logBks;
    }
    //printf("ai done : %lld, in[i] : %lld, approx : %lld\n", ai, in[i], approx);
  }
}
