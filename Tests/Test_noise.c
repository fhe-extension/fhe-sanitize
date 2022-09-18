//#pragma once
#include <stdio.h>	

#include "../random.h"


void myprint(double number)
{
#ifdef _WIN32
	printf("%f\n", number);
#else
	printf("%f\n", number);
#endif
}

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

	double input;
	printf("Choose param : ");
	scanf("%lf", &input);

	long double param = (long double) input;

	printf("Single :\n");

	long double random;
	noise(&random, param);

	myprint((double)random);

	printf("Vector :\n");
	long double *randomv = (long double *) malloc(numSamples * sizeof(long double));
	noise_vector(randomv, param, numSamples);
	
	size_t i;
	long double count = 0;
	for (i = 0; i < numSamples; ++i)
	{
		myprint((double)randomv[i]);
		count += randomv[i];
	}

	free(randomv);
		
	printf("Done !\n");

	totalup((double)count,numSamples);

	printf("Now, odd number of samples\n");

	++numSamples;

	randomv = (long double *) malloc(numSamples * sizeof(long double));
	noise_vector(randomv, param, numSamples);
	
	for (i = 0; i < numSamples; ++i)
	{
		myprint((double)randomv[i]);
		count += randomv[i];
	}

	free(randomv);
		
	printf("Done !\n");

	totalup((double)count,numSamples);

	#ifdef _WIN32
	system("pause"); // Pauses to actually see the output in windows
	#endif

	return 0;
}


