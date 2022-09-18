#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>

#include "../clock.h"
#include "../tgsw.h"
//#include "../tlwe.h"
//#include "../lwe.h"
//#include "fourier.h"

//#include "Test_parameters.h"
//#define TRIES 1e2
#define TEST_EXT_PROD 1
#define TEST_BLINDROTATE 1

int main(int argc, char const *argv[])
{
	srand(time(NULL));

	const int MSG_SPACE = 4;
	const int _lwe_size1 = 100;

	if (TEST_EXT_PROD)
	{
		tlwe_sk sk = tlwe_keygen();

		tgsw_sample tgsw;
		tgsw_sample_init(&tgsw);

		tlwe_sample tlwe;
		tlwe_sample_init(&tlwe);



		int a = rand()%MSG_SPACE, b = rand()%MSG_SPACE;
		//int a = 1, b = 1;


		printf("a = %u ; b = %u\n\n", a, b);
		//system("pause");

		int p[_N];
		size_t i;
		for (i = 1; i < _N; ++i)
			p[i] = 0;


		p[0] = a;
		tgsw_encrypt_over(tgsw, sk, 100.0, p);
		p[0] = b;
		tlwe_encrypt_over(tlwe, sk, 100.0, MSG_SPACE*MSG_SPACE, p);

		tgsw_decrypt_over_and_keep(p, sk, tgsw);
		assert(p[0] == a);

/*		for (size_t r =0; r < _ell * (tlwe_k+1); ++r)
		{
			for (size_t c =0; c < tlwe_k+1; ++c)
			{
				for (size_t x =0; x < tlwe_N; ++x)
					printf("%I64u ", tgsw[x + tlwe_N * (c + r * (tlwe_k + 1))]);
				printf("   ");
			}
			printf("\n");
		}
		system("pause");*/

		tlwe_decrypt_over_and_keep(p, sk, MSG_SPACE*MSG_SPACE, tlwe);
		assert(p[0] == b);

		tlwe_sample tlwe_out = external_product(tlwe, tgsw, 10*_Bg);

		tlwe_decrypt_over_and_keep(p, sk, MSG_SPACE*MSG_SPACE, tlwe_out);

		printf("%u == %u\n", p[0], a*b);

		tgsw_sample_clear(tgsw);
		tlwe_sample_clear(tlwe);
	}

	if (TEST_BLINDROTATE)
	{			
		lwe_sk ns = lwe_keygen(_lwe_size1);
/*		for (int i = 0; i < _lwe_size1; ++i)
			ns[i] = 0;*/

		tlwe_sk S = tlwe_keygen();

		lwe_sk Ns = tlwe_key_extract(S);

		

		int m,m2;
		m = rand()% 2;
		//size_t i;
		ksk ksk = generate_ksk(Ns, ns, 100.0, _k * _N, _lwe_size1);
	    
		/*ksk = (lwe_sample*) malloc(tlwe_k * tlwe_N * PARAM_t * sizeof(lwe_sample));

		for(i = 0; i < tlwe_k * tlwe_N * PARAM_t; ++i) {
			lwe_sample_init(&ksk[i], _lwe_size1);
		}*/

		/*tgsw_sample* bsk;
		bsk = (tgsw_sample*) malloc( _lwe_size1* sizeof(tgsw_sample));
		for(i = 0; i < _lwe_size1; ++i) {
			tgsw_sample_init(bsk+i,tlwe_k,tlwe_N,_ell);
		}*/

		

		bsk bsk = generate_bsk(ns, S, 100.0, _lwe_size1);

/*		int ms[1024];

		for (int i = 0; i < _lwe_size1; ++i)
		{
			tgsw_decrypt(ms, S, bsk[i], tlwe_k, tlwe_N);
			printf("%u == %u\n", ms[0], ns[i]);
		}
*/

		lwe_sample ct = lwe_encrypt(ns, MSG_SPACE, 100.0, m, _lwe_size1);
	    
		bootstrap(ct, bsk, ksk, 10*_Bg, MSG_SPACE,_lwe_size1);


		m2 = lwe_decrypt(ns, MSG_SPACE, ct,_lwe_size1);

		

		printf("%u == %u\n", m, m2);

			//printf("Time for generating ksk:    %f microseconds\n", (double) accum_ksk_gen);// / TRIES);
			//printf("Time for generating ksk:    %f microseconds\n", (double) accum_bs);

		
		lwe_sk_clear(ns);
		lwe_sk_clear(Ns);
		tlwe_sk_clear(S);
		
		ksk_clear(ksk, _k * _N, _lwe_size1);
		bsk_clear(bsk, _lwe_size1);
	}



	#ifdef _WIN32
	system("pause"); // Pauses to actually see the output in windows
	#endif
	
	return 0;
}