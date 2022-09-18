#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

#include "../tgsw.h"
//#include "../tlwe.h"
#include "../clock.h"
//#include "Test_parameters.h"

#define TRIES 1e4


#define N_TRIES 10

#define TLWE_MSG_SPACE 10000

#define TEST_ENC_DEC 0


int main(int argc, char const *argv[])
{
	// allocate memory and create FFT 'plans'
	double accum_rnd = 0, accum_enc = 0, accum_dec = 0;
	long double param = (uint64_t) 1 << (64 - _logBg - 4);

	tlwe_sk sk = tlwe_keygen();

	tgsw_sample ct;
	tgsw_sample_init(&ct);

	int* m= (int*)malloc(_N * sizeof(int));

	for(int j = 0; j < _N; ++j)
			 m[j] = rand() % _Bg;
	
	int* mm= (int*) malloc(_N * sizeof(int));
	for(int try = 0; try < N_TRIES; ++try) {
		start_chrono();

		accum_rnd += stop_chrono();
		
		start_chrono();
		tgsw_encrypt_over(ct, sk, param, m);
		


		accum_enc += stop_chrono();

		start_chrono();
		tgsw_decrypt_over_and_keep(mm, sk, ct);


		accum_dec += stop_chrono();

		// compare the polynomials	
		for (int i = 0; i < _N; ++i) {
			assert(m[i] == mm[i]);
		}
	}

	printf("Time for generating a random polynomial: %f microseconds\n", (double) accum_rnd / N_TRIES);
	printf("Time for encryption:                     %f microseconds\n", (double) accum_enc / N_TRIES);
	printf("Time for decryption:                     %f microseconds\n", (double) accum_dec / N_TRIES);

	// cleanup
	free(m);
	free(mm);
	tlwe_sk_clear(sk);
	tgsw_sample_clear(ct);
	
	#ifdef _WIN32
	system("pause"); // Pauses to actually see the output in windows
	#endif
	
	return 0;
}
