#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "../clock.h"
#include "../fft.h"
//#include "../tlwe.h"
//#include "../lwe.h"
//#include "../random.h"
//#include "Test_parameters.h"

#define N_TRIES 1

#define TLWE_MSG_SPACE 2



int main(int argc, char const *argv[])
{
	fft_init();
	//double accum_enc = 0, accum_dec = 0, accum_ext = 0;
	double accum_sanitize = 0;
	tlwe_sk_fft sk = tlwe_keygen_fft();
	//pkg PK = sanitize_pk_gen_fft(sk, 0.0);

	int* mu= (int*)malloc(_N * sizeof(int));
	
	size_t j;
for (size_t count = 0; count < N_TRIES; ++count){
	for(j = 0; j < _N; ++j)
		mu[j]=  rand()% TLWE_MSG_SPACE;

	tlwe_sample ct = tlwe_encrypt_fft(sk, 0, TLWE_MSG_SPACE, mu);

	//sanitize_pkg_init(P,m,tlwe_k, tlwe_N);


			//epk UPK;
			//sanitize_epk_init(&UPK,tlwe_k,tlwe_N);
		/*	UPK=sanitize_epk_encrypt(PK,m,0,tlwe_k, tlwe_N);
			
			for(j = 0; j < tlwe_N*(tlwe_k+1); ++j)
				ct[j]+=UPK[j];*/

	start_chrono();

	//sanitize_pk_enc(ct, PK);
	sanitize_pk_enc_online_fft(ct, sk, powl(2.0L,22.0L));

	accum_sanitize += stop_chrono();
	int* mmu = tlwe_decrypt_fft(sk, TLWE_MSG_SPACE, ct);
	
			
	for (j = 0; j < _N; ++j)
		assert(mu[j]%TLWE_MSG_SPACE == mmu[j]%TLWE_MSG_SPACE);
		//printf("%u == %u\n", mu[j],mmu[j]);
	printf("\n");
}		
	
	printf("Time for adding a public key encryption of 0 to a TLWE sample: %f microseconds\n", (double) accum_sanitize/N_TRIES);
	
	
	free(mu);

	//sanitize_pk_clear(PK);


	//sanitize_epk_clear(UPK);

	tlwe_sk_clear_fft(sk);
	
	//tlwe_sample_clear(ct);

	#ifdef _WIN32
	system("pause"); // Pauses to actually see the output in windows
	#endif
	
	return 0;




}
