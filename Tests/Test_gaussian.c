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
	printf("The sum is %I64d out of %Iu samples\n",count,numSamples);
#else
	printf("The sum is %lld out of %zu samples\n",count,numSamples);
#endif
}

int main()
{
	size_t numSamples = 1000;

	double input;
	printf("Choose param : ");
	scanf("%lf", &input);

	long double param = (long double) input;
	
	printf("Single :\n");

	uint64_t random;
	gaussian_overZ(&random, param);

	myprint((int64_t)random);

	printf("Vector :\n");
	uint64_t *randomv = (uint64_t *) malloc(numSamples * sizeof(uint64_t));
	gaussian_overZ_vector(randomv, param, numSamples);
	
	size_t i;
	uint64_t count = 0;
	for (i = 0; i < numSamples; ++i)
	{
		myprint((int64_t)randomv[i]);
		count += randomv[i];
	}

	free(randomv);
		
	printf("Done !\n");

	totalup((int64_t)count,numSamples);


	clear_random();

	#ifdef _WIN32
	system("pause"); // Pauses to actually see the output in windows
	#endif

	return 0;
}


