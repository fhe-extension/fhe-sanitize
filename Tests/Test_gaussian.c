//#pragma once
#include <stdio.h>	

#include "../random.h"


void myprint(uint64_t number)
{
	printf("%lld\n", (long long) number);
}

void totalup(uint64_t count, size_t numSamples)
{
	printf("The sum is %" PRId64 " ",count);
        printf("out of %zu samples\n",numSamples);
}

int main()
{
	size_t numSamples = 1000;

	double input;
	printf("Choose param : ");
	scanf("%lf", &input);

	long double param = (long double) input;
        gaussian_param_t p = gaussian(param);

	printf("Single :\n");

	uint64_t random;
	gaussian_overZ(&random, p);

	myprint((int64_t)random);

	printf("Vector :\n");
	uint64_t *randomv = (uint64_t *) malloc(numSamples * sizeof(uint64_t));
	gaussian_overZ_vector(randomv, p, numSamples);
	
	size_t i;
	uint64_t count = 0;
	for (i = 0; i < numSamples; ++i)
	{
		myprint((int64_t)randomv[i]);
		count += randomv[i]*randomv[i];
	}

	free(randomv);
		
	printf("Done !\n");

	totalup((int64_t)count/numSamples,numSamples);


	clear_random();

	#ifdef _WIN32
	system("pause"); // Pauses to actually see the output in windows
	#endif

	return 0;
}


