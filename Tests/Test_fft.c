//#pragma once
#include <stdio.h>
#include <stdlib.h>

#include <complex.h>
#include <fftw3.h>

int main()
{
	size_t N = 1024;
	double *r;
	fftw_complex *c1, *c2;
    fftw_plan p1, p2, b;

    r = (double *) fftw_malloc(sizeof(double) * 2 * N);
    
    c1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N+1));
    c2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N+1));


    if (fftw_import_wisdom_from_filename("../wisdom"))
        {printf("SUCCESS !\n");}

	p1 = fftw_plan_dft_r2c_1d(2 * N, r, c1, FFTW_EXHAUSTIVE | FFTW_DESTROY_INPUT);
    p2 = fftw_plan_dft_r2c_1d(2 * N, r, c2, FFTW_EXHAUSTIVE | FFTW_DESTROY_INPUT);
    b = fftw_plan_dft_c2r_1d(2 * N, c1, r, FFTW_EXHAUSTIVE | FFTW_DESTROY_INPUT);
    fftw_export_wisdom_to_filename("../wisdom");


    size_t i;
    for (i = 0; i < 2 * N; ++i)
    {
    	r[i] = 0;
    }
    r[0] = 3;
    r[1] = 4;

    fftw_execute(p1);

    for (i = 0; i < 2 * N; ++i)
    {
    	r[i] = 0;
    }
    r[0] = 1;
    r[1] = 2;

    fftw_execute(p2);

    for (i = 0; i < N + 1; ++i)
	{
		c1[i] *= c2[i];
	}
	fftw_execute(b);

	for (i = 0; i < 10; ++i)
	{
		printf("%f\n", (double) r[i] / (2 * N));
	}

    fftw_destroy_plan(p1); fftw_destroy_plan(p2); fftw_destroy_plan(b);
    fftw_free(r); fftw_free(c1); fftw_free(c2);

	#ifdef _WIN32
	system("pause"); // Pauses to actually see the output in windows
	#endif

	return 0;
}


