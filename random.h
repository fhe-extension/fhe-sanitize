#pragma once

#include "parameters.h"

#include <stdio.h>
#include <sys/fcntl.h>
//#include <unistd.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>

#ifdef __unix__
#include <fcntl.h>
#include <unistd.h>
#include <sys/random.h>
#endif

#ifdef _WIN32
#include <windows.h>
#endif


//uniform distrib over {0,1}
void random_binary(uint64_t *out);
void random_binary_vector(uint64_t *out, size_t len);

//uniform distrib between 0 and 2^64-1
void uniform64_distribution(uint64_t *out);
void uniform64_distribution_vector(uint64_t *out, size_t len);

//uniform distrib between 0 and 1
void random_double(long double *out);
void random_double_vector(long double *out, size_t len);

//Noise distribution
void noise(long double *out, long double param);
void noise_vector(long double *out, long double param, size_t len);

//Gaussian over Z
void gaussian_overZ(uint64_t *out, long double param);
void gaussian_overZ_vector(uint64_t *out, long double param, size_t len);

//G inverse
void gaussian_ginv(uint64_t *out, uint64_t in, long double param);
void gaussian_ginv_vector(uint64_t *out, uint64_t *in, long double param, size_t len);

//G inverse on polynomials
void gaussian_ginv_poly(uint64_t *out, uint64_t *in, long double param);
void gaussian_ginv_poly_vector(uint64_t *out, uint64_t *in, long double param, size_t len);

//Precomputes stuff
void precompute_random(long double param);
//Clearing precomputation stuff
void clear_random();