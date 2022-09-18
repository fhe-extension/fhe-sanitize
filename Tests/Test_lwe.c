#include <stdio.h>
#include <assert.h>

#include "../clock.h"
//#include "Test_parameters.h"
#include "../lwe.h"


#define n_test 1000
#define TEST_ENC_DEC 1
#define TEST_KS 1

int main(int argc, char const *argv[])
{
	srand(123);

	double accum_enc = 0, accum_dec = 0, accum_gen = 0, accum_ks = 0;

	int m, m_;

	const int _lwe_size1 = 500;
	

	lwe_sample ct;
	lwe_sk skenc = lwe_keygen(_lwe_size1);
	
	lwe_sample_init(&ct, _lwe_size1);
	
/*	for (int i = 0; i < _lwe_size1-1; ++i)
		printf("%I64u\n",skenc[i]);*/
	const int lwe_M = 2;
	const int LWE_SIZE_PRIME = 2048;
	const int LWE_SIZE = 612;
	
	lwe_sk s = lwe_keygen(LWE_SIZE);//cle nouvelle

	lwe_sk ss =	lwe_keygen(LWE_SIZE_PRIME);//cle ancienne

	
	lwe_sample cc;
	lwe_sample_init(&cc, LWE_SIZE_PRIME);//ancien chiffre avec la clé sizeprime

	lwe_sample keyswitched_cc; //nouveau chiffré avec la cle size
	lwe_sample_init(&keyswitched_cc, LWE_SIZE);

	if (TEST_ENC_DEC) 
	{
		for(int i = 0; i < n_test; ++i) 
		{
			//m = (rand()*lwe_M)% PREC -lwe_M/2;
			m = rand()%lwe_M;

			start_chrono();
			lwe_encrypt_over(ct, skenc, lwe_M, 1./(32*lwe_M), m, _lwe_size1);
			//lwe_encrypt(&ct, &skenc, lwe_M, 0, m, _lwe_size1-1);
			accum_enc += stop_chrono();

/*			for (int i = 0; i < _lwe_size1; ++i)
				printf("%I64u\n",ct[i]);*/

			start_chrono();
			m_ = lwe_decrypt_and_keep(skenc, lwe_M, ct, _lwe_size1);
			accum_dec += stop_chrono();

			assert((m%lwe_M) == (m_%lwe_M));
		}
		printf("Time for encryption: %f microseconds\n", (double)accum_enc / n_test);
		printf("Time for decryption: %f microseconds\n", (double)accum_dec / n_test);
	}

	if (TEST_KS) {
		start_chrono();
		ksk ksk = generate_ksk(ss, s, powl(2.0L,49.0L), LWE_SIZE_PRIME, LWE_SIZE);
		accum_gen += stop_chrono();
		for(int try = 0; try < n_test; ++try) {
			m = rand()%lwe_M;
			lwe_encrypt_over(cc, ss, lwe_M, 0.0, m, LWE_SIZE_PRIME);
			
			start_chrono();
			keyswitch_over_and_keep(keyswitched_cc, ksk, LWE_SIZE_PRIME, LWE_SIZE, cc);
			accum_ks += stop_chrono();

			m_ = lwe_decrypt_and_keep(s,lwe_M, keyswitched_cc, LWE_SIZE);
			//assert(m%lwe_M == m_%lwe_M);
			if (m%lwe_M != m_%lwe_M)
				printf("%u, %u\n",m,m_);
		}
		printf("Time for generating ksk: %f microseconds\n", (double)accum_gen / n_test);
		printf("Time for keyswitching:   %f microseconds\n", (double)accum_ks / n_test);
		ksk_clear(ksk, LWE_SIZE_PRIME, LWE_SIZE);
	}
	
	// cleanup
	lwe_sk_clear(skenc);
	lwe_sk_clear(s);
	lwe_sk_clear(ss);
	lwe_sample_clear(ct);
	lwe_sample_clear(cc);
	lwe_sample_clear(keyswitched_cc);

	#ifdef _WIN32
	system("pause"); // Pauses to actually see the output in windows
	#endif
	
	return 0;
}
