//#pragma once
#include <stdio.h>
#include <stdlib.h>

#include <complex.h>
#include <fftw3.h>

int main()
{
	size_t N = 1024;
	long double *r;
	fftwl_complex *c1, *c2;
    fftwl_plan p1, p2, b;

    r = (long double *) fftwl_malloc(sizeof(long double) * 2 * N);
    
    c1 = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * (N+1));
    c2 = (fftwl_complex*) fftwl_malloc(sizeof(fftwl_complex) * (N+1));


    if (fftwl_import_wisdom_from_filename("../wisdom"))
        {printf("SUCCESS !\n");}

	p1 = fftwl_plan_dft_r2c_1d(2 * N, r, c1, FFTW_EXHAUSTIVE | FFTW_DESTROY_INPUT);
    p2 = fftwl_plan_dft_r2c_1d(2 * N, r, c2, FFTW_EXHAUSTIVE | FFTW_DESTROY_INPUT);
    b = fftwl_plan_dft_c2r_1d(2 * N, c1, r, FFTW_EXHAUSTIVE | FFTW_DESTROY_INPUT);
    fftwl_export_wisdom_to_filename("../wisdom");


    size_t i;
    for (i = 0; i < 2 * N; ++i)
    {
    	r[i] = 0;
    }
    r[0] = 3;
    r[1] = 4;

    fftwl_execute(p1);

    for (i = 0; i < 2 * N; ++i)
    {
    	r[i] = 0;
    }
    r[0] = 1;
    r[1] = 2;

    fftwl_execute(p2);

    for (i = 0; i < N + 1; ++i)
	{
		c1[i] *= c2[i];
	}
	fftwl_execute(b);

	for (i = 0; i < 10; ++i)
	{
		printf("%f\n", (double) r[i] / (2 * N));
	}

    fftwl_destroy_plan(p1); fftwl_destroy_plan(p2); fftwl_destroy_plan(b);
    fftwl_free(r); fftwl_free(c1); fftwl_free(c2);

	#ifdef _WIN32
	system("pause"); // Pauses to actually see the output in windows
	#endif

	return 0;
}


