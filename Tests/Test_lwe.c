#include <stdio.h>
#include <assert.h>

#include "../clock.h"
//#include "Test_parameters.h"
#include "../lwe.h"


#define _two64_double powl(2.0L,64.0L)
#define _PI acosl(-1.0L)

#define n_test 100
#define TEST_KS 1

int main(int argc, char const *argv[])
{
	srand(123);
	long double q = powl(2.0L,_ell*_logBg);
	double accum_gen = 0, accum_ks = 0;
    long double accum_var = 0;
	int m, m_;

	long double stddev_ks = powl(2.0L,-13.6L);
	long double stddev_bk = powl(2.0L,-33.8L);
	//lwe_sk skenc = lwe_gaussian_keygen(_lwe_size1, powl(2.0L, 1.2L));
	  
/*	for (int i = 0; i < _lwe_size1-1; ++i)
		printf("%I64u\n",skenc[i]);*/
	const int lwe_M = 2;
	const int LWE_SIZE_PRIME = 1024;
	const int LWE_SIZE = 538;
	
	lwe_sk s = lwe_keygen(LWE_SIZE);// new secret key
	
	lwe_sk ss = lwe_gaussian_keygen(LWE_SIZE_PRIME,sqrtl(2.0L*_PI)*stddev_bk*q);
	
	//lwe_sk ss =lwe_keygen(LWE_SIZE_PRIME);

	lwe_sample cc;
	lwe_sample_init(&cc, LWE_SIZE_PRIME);//old ciphertext with secret key LWE_SIZE_PRIME

	lwe_sample keyswitched_cc; //new ciphertext with secret key LWE_SIZE
	lwe_sample_init(&keyswitched_cc, LWE_SIZE);

	
	if (TEST_KS) {
          accum_var = 0;
          start_chrono();
          //ksk ksk = generate_ksk(ss, s, powl(2.0L,50.4L), LWE_SIZE_PRIME, LWE_SIZE);
          ksk ksk = generate_ksk(ss, s,stddev_ks*_two64_double, LWE_SIZE_PRIME, LWE_SIZE);
          accum_gen += stop_chrono();
          for(int try = 0; try < n_test; ++try) {
            m = rand()%lwe_M;
            //lwe_encrypt_over(cc, ss, lwe_M, (1L<<59)/lwe_M, m, LWE_SIZE_PRIME);
            	lwe_encrypt_over(cc, ss, lwe_M, 0.0L, m, LWE_SIZE_PRIME);
            start_chrono();
            keyswitch_over_and_keep(keyswitched_cc, ksk, LWE_SIZE_PRIME, LWE_SIZE, cc);
            accum_ks += stop_chrono();

            long double err;
            m_ = lwe_decrypt_and_keep(s,lwe_M, keyswitched_cc, LWE_SIZE, &err);
            accum_var += err*err;
            //assert(m%lwe_M == m_%lwe_M);
            if (m%lwe_M != m_%lwe_M)
              printf("%u, %u\n",m,m_);
          }
          printf("Time for generating ksk: %f microseconds\n", (double)accum_gen / n_test);
          printf("Time for keyswitching:   %f microseconds\n", (double)accum_ks / n_test);
          printf("Variance of the error : %f\n", (double)accum_var/n_test);
          ksk_clear(ksk, LWE_SIZE_PRIME, LWE_SIZE);
        }

        // cleanup

        lwe_sk_clear(s);
        lwe_sk_clear(ss);

        lwe_sample_clear(cc);
        lwe_sample_clear(keyswitched_cc);

	#ifdef _WIN32
	system("pause"); // Pauses to actually see the output in windows
	#endif

	return 0;
}
