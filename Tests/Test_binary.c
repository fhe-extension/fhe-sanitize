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

void totalup(uint64_t count, size_t numSamples)
{
#ifdef _WIN32
	printf("I got %Iu ones out of %Iu samples\n",count,numSamples);
#else
	printf("I got %llu ones out of %zu samples\n",count,numSamples);
#endif
}

int main()
{
	size_t numSamples = 1001;

	printf("Single :\n");

	uint64_t random;
	random_binary(&random);

	myprint(random);

	printf("Vector :\n");
	uint64_t *randomv = (uint64_t *) malloc(numSamples * sizeof(uint64_t));
	random_binary_vector(randomv, numSamples);
	
	size_t i;
	int count = 0;
	for (i = 0; i < numSamples; ++i)
	{
		myprint(randomv[i]);
		count += randomv[i];
	}
		
	free(randomv);
	printf("Done !\n");

	totalup(count, numSamples);

	#ifdef _WIN32
	system("pause"); // Pauses to actually see the output in windows
	#endif

	return 0;
}


