//#pragma once
#include <stdio.h>

#include "../random.h"



int main()
{
	size_t numSamples = 1000000;

	double input;
	printf("Choose param : ");
	scanf("%lf", &input);

	double param = (double) input;
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
	int64_t avg = 0;
	for (i = 0; i < numSamples; ++i)
	{
          printf("%lld\n", (long long int) randomv[i]);
          avg += randomv[i];
          count += randomv[i]*randomv[i];
	}

	free(randomv);

	printf("Done !\n");

	printf("Observed average  : %lf\n",(double) avg/numSamples);
	printf("Observed variance : %lf\n",(double) count/numSamples);
        printf("Expected variance : %lf\n", (param*param)/(2*acos(-1.0L)));

        clear_gaussian_param(p);
	clear_random();

	#ifdef _WIN32
	system("pause"); // Pauses to actually see the output in windows
	#endif

	return 0;
}


