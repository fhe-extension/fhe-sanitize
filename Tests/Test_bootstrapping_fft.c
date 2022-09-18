#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>

#include "../clock.h"
#include "../fft.h"

#define _two64_double powl(2.0L,64.0L)

int main(int argc, char const *argv[])
{
	srand(time(NULL));
	fft_init();

	const int MSG_SPACE = 2;
	const int _n = 612;
/*	uint64_t param = 1;
	param <<= 10;*/
	long double param = 30825788;

	precompute_random(param);

	lwe_sk ns = lwe_keygen(_n);
	tlwe_sk_fft S = tlwe_keygen_fft();
	lwe_sk Ns = tlwe_key_extract_fft(S);

	ksk ksk = generate_ksk(Ns, ns, powl(2.0L,49.0L), _k * _N, _n);
	bsk_fft bsk = generate_bsk_fft(ns, S, powl(2.0L,22.0L), _n);
	
	printf("Precompute PkEnc\n");

	pkg PK = sanitize_pk_gen_fft(S, powl(2.0L,22.0L));
	//pkg PK = sanitize_pk_gen_fft(S, 0.0);
	precompute_pkenc(PK);

	int m,m2;
	m = rand()% 2;

	lwe_sample ct = lwe_encrypt(ns, MSG_SPACE, powl(2.0L,57.0L), m, _n);
	   

	printf("GO !\n");

	start_chrono();

	sanitize_fft(ct, bsk, ksk, PK, param, MSG_SPACE,_n);
	//bootstrap_fft(ct, bsk, ksk, param, MSG_SPACE,_n);


	long double temps = stop_chrono();

	printf("STOP !\n");

	m2 = lwe_decrypt(ns, MSG_SPACE, ct,_n);

	printf("%u == %u\n", m, m2);
	
	printf("Time for sanitization: %lf microseconds\n", (double) temps);
		
	ksk_clear(ksk, _k * _N, _n);
	bsk_clear_fft(bsk, _n);


	lwe_sk_clear(ns);
	tlwe_sk_clear_fft(S);
	lwe_sk_clear(Ns);

	sanitize_pk_clear(PK);
	clear_pkenc();
	clear_random();

	fft_clear();

	#ifdef _WIN32
	system("pause"); // Pauses to actually see the output in windows
	#endif
	
	return 0;
}