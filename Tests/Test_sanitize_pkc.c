#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <fcntl.h> // for open
#include <unistd.h> // for close

#include "../clock.h"
#include "../fft.h"
//#include "../lwe.h"
//#include "../random.h"
//#include "Test_parameters.h"

#define N_TRIES 100

#define TLWE_MSG_SPACE 4



int main(int argc, char const *argv[])
{
	double accum_sanitize = 0;
	srand(time(NULL));


	fft_init();

	tlwe_sk_fft tsk = tlwe_keygen_fft(powl(2.0L,1.2L));
	tlwe_sample ct;
	tlwe_sample_init(&ct);

	pkc_fft PK = sanitize_pkc_gen_fft(tsk, sqrtl(2.0L*acosl(-1.0L))*powl(2.0L,1.2L));


	printf("PK generated\n");

	#ifdef _WIN32
	system("pause"); // Pauses to actually see the output in windows
	#endif


	int* mu= (int*)malloc(_N * sizeof(int));
	int* mmu= (int*)malloc(_N * sizeof(int));
	size_t j;
	for(j = 0; j < _N; ++j)
		mu[j]=  rand()% TLWE_MSG_SPACE;

	tlwe_encrypt_over_fft(ct, tsk, 0, TLWE_MSG_SPACE, mu);

	//tlwe_sample_zero(ct);

	printf("Sanitizing\n");

	start_chrono();

	sanitize_pkc_enc_fft(ct, PK, gaussian(sqrtl(2.0L*acosl(-1.0L))*powl(2.0L,8.4L)),gaussian(sqrtl(2.0L*acosl(-1.0L))*powl(2.0L,8.4L)),gaussian(sqrtl(2.0L*acosl(-1.0L))*powl(2.0L,15.95L)));

	accum_sanitize += stop_chrono();
	printf("Sanitized\n");
	tlwe_decrypt_over_and_keep_fft(mmu, tsk, TLWE_MSG_SPACE, ct);
	
			
	for (j = 0; j < _N; ++j)
          if (mu[j] % TLWE_MSG_SPACE != mmu[j] % TLWE_MSG_SPACE)
            printf("%u == %u\n",mu[j],mmu[j]);
        printf("\n");
//}		
	
	printf("Time for adding a public key encryption of 0 to a TLWE sample: %f microseconds\n", (double) accum_sanitize/N_TRIES);
	
	free(mu);
	free(mmu);
	//sanitize_pkc_clear_fft(PK);
	
	tlwe_sk_clear_fft(tsk);
	
	tlwe_sample_clear(ct);

	fft_clear();

	#ifdef _WIN32
	system("pause"); // Pauses to actually see the output in windows
	#endif
	
	return 0;




}
