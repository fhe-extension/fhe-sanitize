#include "clock.h"
#ifdef _WIN32
struct timespec start, stop;
#else
struct timeval start, stop;
#endif

void start_chrono()
{
#ifdef _WIN32
	clock_gettime(CLOCK_REALTIME, &start);
#else
	gettimeofday(&start, NULL);
#endif
}

double stop_chrono()
{
#ifdef _WIN32
	clock_gettime(CLOCK_REALTIME, &stop);
	return ((stop.tv_sec*1e9 + stop.tv_nsec) - (start.tv_sec*1e9 + start.tv_nsec))/1000.0;
#else
	gettimeofday(&stop, NULL);
	return (stop.tv_sec*1e6 + stop.tv_usec) - (start.tv_sec*1e6 + start.tv_usec);
#endif
}