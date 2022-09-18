//#pragma once
#include <stdio.h>	

#include "../random.h"


void myprint(uint64_t number)
{
#ifdef _WIN32
	printf("%I64d\n", number);
#else
	printf("%lld\n", number);
#endif
}

void totalup(uint64_t count, size_t numSamples)
{
#ifdef _WIN32
	printf("The sum is %I64u out of %Iu samples\n",count,numSamples);
#else
	printf("The sum is %llu out of %zu samples\n",count,numSamples);
#endif
}

int main()
{
	uint64_t n = 3;

	uint64_t *v = (uint64_t *) malloc(n*sizeof(uint64_t));

	v[0] = (uint64_t) 17 << (64 - _logBg * _ell);
	v[1] = (uint64_t) 54 << (64 - _logBg * _ell);
	v[2] = (uint64_t) 107 << (64 - _logBg * _ell);

	double input;
	printf("Choose param : ");
	scanf("%lf", &input);

	long double param = (long double) input;

	uint64_t *random = (uint64_t *) malloc(_ell*n*sizeof(uint64_t));
	gaussian_ginv_vector(random, v, param, n);

/*	uint64_t u = 1;
	uint64_t count = 0;
	size_t i;
	size_t j;
	for (j=0; j<n; j++)
	{
		for (i = ell-1; i + 1 > 0; --i)
		{
			myprint((int64_t)random[j*ell+i]);
			count += random[j*ell+i]*u;
			u*=Bg;
		}

	totalup(count,ell);
	u = 1;
	count = 0;
	}*/
	
	uint64_t tmp = (uint64_t) 1 << (64 - _logBg);
	uint64_t count = 0;
	for (size_t i = 0; i < n; ++i)
	{
		for (size_t power = 0; power < _ell; ++power)
		{
			count += random[i*_ell+power] * tmp;
			tmp >>= _logBg;
		}

		totalup(count >> (64 - _logBg * _ell), _ell);
		tmp = (uint64_t) 1 << (64 - _logBg);
		count = 0;
	}
	
	clear_random();

	#ifdef _WIN32
	system("pause"); // Pauses to actually see the output in windows
	#endif


	return 0;
}


