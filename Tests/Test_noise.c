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
	size_t numSamples = 1000000;

	double input;
	printf("Choose param : ");
	scanf("%lf", &input);

	double param = (double) input;

	printf("Single :\n");

	double random;
	noise(&random, param);

	myprint((double)random);

	printf("Vector :\n");
	double *randomv = (double *) malloc(numSamples * sizeof(double));
	noise_vector(randomv, param, numSamples);
	
	size_t i;
	double count = 0;
	for (i = 0; i < numSamples; ++i)
	{
		myprint((double)randomv[i]);
		count += roundl(randomv[i])*roundl(randomv[i]);
	}

	free(randomv);
		
	printf("Done !\n");

        double real_variance = (double)count/numSamples;
	totalup(real_variance,numSamples);
        printf("Expected : %lf, Difference : %lf, Relative Difference %lf\n", (double) (param*param),(double) (real_variance-param*param),(double) ((real_variance-param*param)/real_variance));

	printf("Now, odd number of samples\n");

	++numSamples;

	randomv = (double *) malloc(numSamples * sizeof(double));
	noise_vector(randomv, param, numSamples);
	
	for (i = 0; i < numSamples; ++i)
	{
		//myprint((double)randomv[i]);
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


