//#pragma once
#include <stdio.h>	

#include "../random.h"


void myprint(uint64_t number)
{
#ifdef _WIN32
	printf("%Iu\n", number);
#else
	printf("%llu\n", number);
#endif
}

int main()
{
	size_t numSamples = 1001;

	printf("Single :\n");

	uint64_t random;

	uniform64_distribution(&random);

	myprint(random);

	printf("Vector :\n");
	uint64_t *randomv = (uint64_t *) malloc(numSamples * sizeof(uint64_t));
	uniform64_distribution_vector(randomv, numSamples);
	
	size_t i;
	for (i = 0; i < numSamples; ++i)
	{	
		myprint(randomv[i]);
	}
	free(randomv);
	printf("Done !\n");

	#ifdef _WIN32
	system("pause"); // Pauses to actually see the output in windows
	#endif

	return 0;
}


