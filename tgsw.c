#include "tgsw.h"

#define _two64_double powl(2.0L,64.0L) 

void tgsw_sample_init(tgsw_sample *ct)
{
	*ct = (tgsw_sample) malloc((_k+1) * (_k+1) * _N * _ell * sizeof(uint64_t));
}

void tgsw_sample_clear(tgsw_sample ct)
{
	free(ct);
}

void tgsw_sample_zero(tgsw_sample *ct)
{
	size_t max = (_k+1)*(_k+1)*_N*_ell;
	for(size_t i=0; i < max; ++i) {
		(*ct)[i] = 0;
	}
}

/* G SPECS */
/* (d+1)ell ROWS, (d+1)N COLUMNS, Each adjacent N numbers are polynomial coefficients, Powers of Bg decreasing top to bottom from Bg^(ell-1) to 1 */
//_ell and _Bg hard coded
tgsw_sample tgsw_encrypt(tlwe_sk sk, long double param, int *m)
{
	tgsw_sample ct;
	tgsw_sample_init(&ct);

	size_t i;
	for (i = 0; i < (_k+1) * _ell; ++i)
		tlwe_encrypt_zero_over(ct + i * ( _k + 1 ) * _N, sk, param);

	uint64_t tmp = (uint64_t) 1 << (64 - _logBg);

	size_t coeff;
	for (size_t power = 0; power < _ell; ++power)
	{
		for (i = 0; i <= _k; ++i)
		{
			for (coeff = 0; coeff < _N; coeff++)
			ct[power * (_k+1) * _N + i * (_ell * (_k+1) + 1) * _N + coeff] += m[coeff] * tmp;
		}
		tmp >>= _logBg;
	}
	return ct;
}
void tgsw_encrypt_over(tgsw_sample ct, tlwe_sk sk, long double param, int *m)
{
	size_t i;
	for (i = 0; i < (_k+1) * _ell; ++i)
		tlwe_encrypt_zero_over(ct + i * ( _k + 1 ) * _N, sk, param);

	uint64_t tmp = (uint64_t) 1 << (64 - _logBg);

	size_t coeff;
	for (size_t power = 0; power < _ell; ++power)
	{
		for (i = 0; i <= _k; ++i)
		{
			for (coeff = 0; coeff < _N; coeff++)
			ct[power * (_k+1) * _N + i * (_ell * (_k+1) + 1) * _N + coeff] += m[coeff] * tmp;
		}
		tmp >>= _logBg;
	}
}


int* tgsw_decrypt(tlwe_sk sk, tgsw_sample ct)
{
	int *m = (int *) malloc(_N * sizeof(int));
	tlwe_decrypt_over_and_keep(m, sk, _Bg, ct + (_k+1) * _N * _ell * _k);
	tgsw_sample_clear(ct);
	return m;
}
void tgsw_decrypt_over_and_keep(int *m, tlwe_sk sk, tgsw_sample ct)
{
	tlwe_decrypt_over_and_keep(m, sk, _Bg, ct + (_k+1) * _N * _ell * _k);
}

bsk generate_bsk(lwe_sk sk_in, tlwe_sk sk_out, long double param, size_t n)
{
	bsk BSK = (tgsw_sample *) malloc(n * sizeof(tgsw_sample));;
	int* p = (int *) malloc(_N * sizeof(int));
	size_t i;
	for (i = 1; i < _N; ++i)
		p[i] = 0;

	for (i = 0; i < n; ++i)
	{
		p[0] = sk_in[i];
		BSK[i] = tgsw_encrypt(sk_out, param, p);
	}
	free(p);

	return BSK;
}

void bsk_clear(bsk bsk, size_t n)
{
	for (size_t i = 0; i < n; ++i)
		tgsw_sample_clear(bsk[i]);
	free(bsk);
}

tlwe_sample external_product(tlwe_sample tlwe, tgsw_sample tgsw, long double param)
{
	tlwe_sample out;
	tlwe_sample_init(&out);
	uint64_t* decomposed = (uint64_t *) malloc((_k+1) * _ell * _N * sizeof(uint64_t));
	size_t i, j;

	gaussian_ginv_poly_vector(decomposed, tlwe, param, _k+1);

	for(i = 0; i < (_k+1)*_N; ++i) 
			out[i] = 0;
	for(i = 0; i < (_k+1); ++i) // Column number in the matrix
		for (j = 0; j<(_k+1)*_ell; ++j) // Row number in the matrix
			multiply_accumulate_poly(out + i * _N, decomposed + j * _N, tgsw + (i + j * (_k+1)) * _N);

	free(decomposed);
	return out;
}

void accum_external_product(tlwe_sample out, tlwe_sample tlwe, tgsw_sample tgsw, long double param)
{
	uint64_t* decomposed = (uint64_t *) malloc((_k+1) * _ell * _N * sizeof(uint64_t));
	size_t i, j;

	gaussian_ginv_poly_vector(decomposed, tlwe, param, _k+1);

	for(i = 0; i < (_k+1); ++i) // Column number in the matrix
		for (j = 0; j<(_k+1)*_ell; ++j) // Row number in the matrix
			multiply_accumulate_poly(out + i * _N, decomposed + j * _N, tgsw + (i + j * (_k+1)) * _N);

	free(decomposed);
}


tlwe_sample blind_rotate(lwe_sample lwe, long double param, int M, size_t n, uint64_t* testv, tgsw_sample *bsk)
{		
	size_t max=_k*_N; 
	size_t _2N=2*_N;
	
	size_t i, j, pow;	

	tlwe_sample out, tmp_acc;
	tlwe_sample_init(&out);
	tlwe_sample_init(&tmp_acc);

	/*int * ct_bar = (int *) malloc((n+1) * sizeof(int));
	//size_t i, pow;//, j,pow;
	//copie de lwe[i] dans ct_bar
	for( i = 0; i <= n; ++i) {
		ct_bar[i] = (int) (lwe[i] >> (64 - _log2N));
			printf("%d\n",ct_bar[i]);
	}
	int bbar=ct_bar[n];*/

/*	uint64_t* testvbis = (uint64_t*) malloc(N * sizeof(uint64_t));
	for (i = 0; _2N + i - bbar < N; ++i)
		testvbis[i] = testv[_2N+i-bbar];
	for (; N + i - bbar < N; ++i)
		testvbis[i] = -testv[N+i-bbar];
	for (; i < N; ++i)
		testvbis[i] = testv[i-bbar];*/

	int bbar = (lwe[n] >> (64 - _log2N));
	for (i = 0; i < max; ++i)
		out[i] = 0;

	if (bbar < _N)
	{
		for (i = 0; i < bbar; ++i)
			out[max+i] = -testv[_N+i-bbar];
		for (; i < _N; ++i)
			out[max+i] = testv[i-bbar];
	}
	else
	{
		for (i = 0; i < bbar - _N; ++i)
			out[max+i] = testv[_2N-bbar+i];
		for (; i < _N; ++i)
			out[max+i] = -testv[_N-bbar+i];
	}
/*	for (i = 0; _2N + i - bbar < N; ++i)
		out[max+i] = testv[_2N+i-bbar];
	for (; N + i - bbar < N; ++i)
		out[max+i] = -testv[N+i-bbar];
	for (; i < N; ++i)
		out[max+i] = testv[i-bbar];*/

/*	printf("bbar = %u\n", bbar);
	for (i = 0; i < N; ++i)
		printf("%I64d ", out[max+i]);*/

/*	int bbar=ct_bar[n_in];

	printf("%d\n", bbar);
	//construction de Xbar
	uint64_t* Xbbar =(uint64_t*) malloc(N * sizeof(uint64_t));
	
	for( i = 0; i < N; ++i) 
		Xbbar[i] = 0;

	if(bbar < N)
		Xbbar[bbar] = 1;
	else
		Xbbar[bbar - N] = -1;
	
	printf("BLAH5 \n");
	system("pause");


	//testv=testv*X^bbar
	uint64_t* testvbis=(uint64_t*) malloc(N * sizeof(uint64_t));

	for(i = 0; i < N; ++i) {
		testvbis[i] = 0;
	}
	multiply_accumulate_poly(testvbis,testv,Xbbar,N);
	printf("BLAH7 \n");
	memcpy(out+max,testvbis,N);


	memcpy(out+max,testvbis,N);*/

/*	tlwe_sample tmp_acc; 
	tlwe_sample_init(&tmp_acc,k,N);
	tlwe_sample_zero(tmp_acc,k,N);
	//initialisation de acc
	memcpy(out+max, testvbis, N);*/
	

/*	for (i = 0; i < n_in;  ++i) {
		if(ct_bar[i]!=0) {	
			printf("BLAH8 \n");
			
			external_product(tmp_acc,out,bsk[i], param, k, N);
		
			for(pow = 0; pow < N; ++pow) {
				testvbis[pow] = 0;
			}
			testvbis[0] = -1;
			
			if (ct_bar[i] < N)
				testvbis[N-ct_bar[i]] -= 1;
			else
				testvbis[2*N-ct_bar[i]] += 1;

			tlwe_add_mul(out,tmp_acc,testvbis,k,N);
        
		}
	}*/

	int abari;
	for (i = 0; i < n;  ++i)
	{

//		printf("BLAHBLAH %u\n", i);
		abari = (-lwe[i] >> (64 - _log2N));
//		printf("abari : %u\n", abari);
//		system("pause");
		if (abari < _N)
		{
			for (pow = 0; pow < abari; ++pow)
				for (j = 0; j <= _k; ++j)
					tmp_acc[j * _N + pow] = - out[j * _N + pow] - out[j * _N + _N + pow - abari];
			for (; pow < _N; ++pow)
				for (j = 0; j <= _k; ++j)
					tmp_acc[j * _N + pow] = - out[j * _N + pow] + out[j * _N + pow - abari];
		}
		else
		{
			for (pow = 0; pow < abari - _N; ++pow)
				for (j = 0; j <= _k; ++j)
					tmp_acc[j * _N + pow] = - out[j * _N + pow] + out[j * _N + _2N + pow - abari];
			for (; pow < _N; ++pow)
				for (j = 0; j <= _k; ++j)
					tmp_acc[j * _N + pow] = - out[j * _N + pow] - out[j * _N + _N + pow - abari];
		}

		/*for (pow = 0; _2N + pow - abari < N; ++pow)
			for (j = 0; j <= k; ++j)
				tmp_acc[j * N + pow] = - out[j * N + pow] + out[j * N + _2N + pow - abari];
		for (; N + pow - abari < N; ++pow)
			for (j = 0; j <= k; ++j)
				tmp_acc[j * N + pow] = - out[j * N + pow] - out[j * N + N + pow - abari];
		for (; pow < N; ++pow)
			for (j = 0; j <= k; ++j)
				tmp_acc[j * N + pow] = - out[j * N + pow] + out[j * N + pow - abari];*/

//		printf("Before ext prod\n");
//		system("pause");

		accum_external_product(out, tmp_acc, bsk[i], param);
	}


	//free(ct_bar);
	//free(Xbbar);
	//free(testvbis);
	//lwe_sample_clear(lwe);
	tlwe_sample_clear(tmp_acc);

	return out;
}



void bootstrap(lwe_sample ct, bsk bsk, ksk ksk, long double param, int M_out, size_t n)
{
	uint64_t* test_vector=(uint64_t *) malloc(_N * sizeof(uint64_t));
   
	size_t i;
	for (i = 0; i < _N; ++i)
		test_vector[i] = - _two64_double / (2 * M_out);

/*	size_t bound = N/2;
	size_t i;
	for (i = 0; i < bound; ++i)
		test_vector[i] = - _two64_double / (2 * M_out);

		// test_vector.coeffs[i] = ((int) ((double)(i*msg_space_in)/(2*TLWE_N)+0.5)) / (double) msg_space_out;
		//test_vector[i] = (-1*(_two64_double)) / M_out;
		// printf("[%4d]  %f\n", i, test_vector.coeffs[i]);
	for (; i < N; ++i)
		test_vector[i] = _two64_double / (2 * M_out);
	*/
		
/*	for (i = 0; i < halfN; ++i)
		test_vector[i] = - _two64_double/M_out;
	for (; i < N; ++i)
		test_vector[i] = _two64_double/M_out;*/
	printf("Entering blind_rotate...\n");
	
	tlwe_sample tmp_tlwe = blind_rotate(ct, param, M_out,n,test_vector,bsk);

	tmp_tlwe[_k * _N] += _two64_double / 8;

	// sample extract
	printf("Entering tlwe_extract...\n");
	
	lwe_sample xtracted_lwe = tlwe_extract(tmp_tlwe);

	printf("Entering keyswitch...\n");
	// key switch
	keyswitch_over_and_keep(ct, ksk,_k*_N,n,xtracted_lwe);

	tlwe_sample_clear(xtracted_lwe);

	//out[n] += _two64_double / 8;

	// cleanup
	free(test_vector);
}
