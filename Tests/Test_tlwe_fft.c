#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "../clock.h"
#include "../fft.h"
//#include "../lwe.h"
//#include "../random.h"

#define N_TRIES 10000

#define TLWE_MSG_SPACE 4

#define TEST_ENC_DEC 1
#define TEST_EXTRACT 1


int main(int argc, char const *argv[])
{
	fft_init();
	//double accum_rnd = 0, 
	double accum_enc = 0, accum_dec = 0, accum_ext = 0;

	tlwe_sk_fft sk = tlwe_keygen_fft(powl(2.0L,1.2L));

	int* m= (int*)malloc(_N * sizeof(int));
	int* mm= (int*) malloc(_N * sizeof(int));
	//poly_init(&m,tlwe_N);


	tlwe_sample ct;
	tlwe_sample_init(&ct);
	//tlwe_sample_zero(&ct,tlwe_k,tlwe_N);

	

	//For extraction
	lwe_sample ext_ct;

	lwe_sample_init(&ext_ct,_N*_k);

	lwe_sk ext_key = tlwe_key_extract_fft(sk);
	int ext_m;
	

	if (TEST_ENC_DEC) {
	//	for (int try = 0; try < N_TRIES; ++try) {
			
			start_chrono();
	//		int N=10; int k=2; 
	//		int max=N*k;
	        
	      
			for(int j = 0; j < _N; ++j)
			 {
				m[j]=  rand()% TLWE_MSG_SPACE;
				printf("%d", m[j]);
				printf(" "); 
			 }
			printf("\n");
			//print_polynomial(m,tlwe_N);
			

            //tlwe_encrypt(ct, sk, 1./(16*TLWE_MSG_SPACE), TLWE_MSG_SPACE, m,tlwe_k,tlwe_N);
			tlwe_encrypt_over_fft(ct, sk, powl(2.0L,57.0L), TLWE_MSG_SPACE, m);
			accum_enc += stop_chrono();

			
			double *err_N = malloc(_N*sizeof(double));
			//gettimeofday(&_begin, NULL);
			start_chrono();
			tlwe_decrypt_over_and_keep_fft(mm, sk, TLWE_MSG_SPACE, ct, err_N);
			for(int i = 0; i < _N; ++i) 
			{
       			printf("%d", mm[i]);
				printf(" "); 
			}
			printf("\n"); printf("\n");printf("\n"); 
			accum_dec += stop_chrono();
                        free(err_N);

			for (size_t i = 0; i < _N; ++i) 
			{
				assert((m[i]%TLWE_MSG_SPACE) == (mm[i]%TLWE_MSG_SPACE));
			}
			printf("\n"); 
		}
		//printf("Time for generating a random polynomial: %f microseconds\n", (double) accum_rnd / N_TRIES);
		//printf("Time for encryption: %f microseconds\n", (double) accum_enc / N_TRIES);
		//printf("Time for decryption: %f microseconds\n", (double) accum_dec / N_TRIES);
	if (TEST_EXTRACT) {
          double stddev = powl(2.0L, -14.0L);
          double var = 0.;
          for (int try = 0; try < N_TRIES; ++try) {
            srand((try%245)*11);


            for (int j = 0; j < _N; ++j) {
              m[j] = rand()% TLWE_MSG_SPACE;
            }
            

            tlwe_encrypt_over_fft(ct, sk, stddev * powl(2.0L,64.0L), TLWE_MSG_SPACE, m);
			
            start_chrono();
			

			
            tlwe_extract_over_and_keep(ext_ct, ct);
            accum_ext += stop_chrono();
			
            //n --> n-1=k*N
            double err;
            ext_m = lwe_decrypt_and_keep(ext_key, TLWE_MSG_SPACE, ext_ct, _k * _N, &err);
            var += err * err;
            //int rounded= m[0];
            //printf("%d", rounded); printf(" and "); printf("%d", ext_m); printf("\n");
			
            assert((m[0]%TLWE_MSG_SPACE)==(ext_m%TLWE_MSG_SPACE));
          }
          printf("Time for extracting an LWE sample: %f microseconds\n", (double) accum_ext/N_TRIES);
          printf("Measured variance is : %lf\n", var / N_TRIES);
          printf("Expected variance is : %lf\n", stddev * stddev);
          printf("Relatice difference : %lf\n", var / N_TRIES / (stddev * stddev) - 1);
	}

	
	free(m);
	free(mm);
	tlwe_sk_clear_fft(sk);
	lwe_sk_clear(ext_key);
	tlwe_sample_clear(ct);
        lwe_sample_clear(ext_ct);

	#ifdef _WIN32
	system("pause"); // Pauses to actually see the output in windows
	#endif
	
	return 0;




}
