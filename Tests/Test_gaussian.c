//#pragma once
#include <stdio.h>

#include "../random.h"



int main()
{
	size_t numSamples = 1000000;

	double input;
	printf("Choose param : ");
	scanf("%lf", &input);

	long double param = (long double) input;
        gaussian_param_t p = gaussian(param);

	printf("Single :\n");

	uint64_t random;
	gaussian_overZ(&random, p);

        printf("%llu\n", (long long int) random);

	printf("Vector :\n");
	uint64_t *randomv = (uint64_t *) malloc(numSamples * sizeof(uint64_t));
	gaussian_overZ_vector(randomv, p, numSamples);

	size_t i;
	uint64_t count = 0;
	for (i = 0; i < numSamples; ++i)
	{
		printf("%llu\n", (long long int) randomv[i]);
		count += randomv[i]*randomv[i];
	}

	free(randomv);

	printf("Done !\n");

	printf("Observed variance : %Lf\n",(long double) count/numSamples);
        printf("Expected : %Lf\n", (param*param)/(2*acos(-1.0L)));

	clear_random();

	#ifdef _WIN32
	system("pause"); // Pauses to actually see the output in windows
	#endif

	return 0;
}


