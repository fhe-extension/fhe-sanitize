//#pragma once
#include <stdio.h>	

#include "../random.h"

// NOTE : I cast the random elements to double because I was unable to easily print double values.
void totalup(double count, size_t numSamples)
{
#ifdef _WIN32
	printf("The sum is %f out of %Iu samples\n",count,numSamples);
#else
	printf("The sum is %f out of %zu samples\n",count,numSamples);
#endif
}

int main()
{
	size_t numSamples = 1000;

	printf("Single :\n");

	double random;
	random_double(&random);

	printf("%f\n", (double) random);

	printf("Vector :\n");
	double *randomv = (double *) malloc(numSamples * sizeof(double));
	random_double_vector(randomv, numSamples);
	
	size_t i;
	double count = 0.;
	for (i = 0; i < numSamples; ++i)
	{
		printf("%f\n",(double) randomv[i]);
		count += randomv[i];
	}
		
	free(randomv);
	printf("Done !\n");

	totalup((double) count,numSamples);

	#ifdef _WIN32
	system("pause"); // Pauses to actually see the output in windows
	#endif

	return 0;
}


