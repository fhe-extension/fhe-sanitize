#pragma once

#include "parameters.h"

#include <stdio.h>
#include <sys/fcntl.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>


#ifdef _WIN32
#include <windows.h>
#else
#include <fcntl.h>
#include <unistd.h>
#include <sys/random.h>
#endif

//Gaussian parameters, including auxilliary inputs to the sampler
struct gaussian_param {
  long double param;
  size_t target;
  uint64_t* z;
  uint64_t adjust;
};
typedef struct gaussian_param gaussian_param_t;

//uniform distrib over {0,1}
void random_binary(uint64_t *out);
void random_binary_vector(uint64_t *out, size_t len);

void random_binary_vector_int64(int64_t *out, size_t len);
void uniform64_distribution_int64(int64_t *out);
void uniform64_distribution_vector_int64(int64_t *out, size_t len);


//uniform distrib between 0 and 2^64-1
void uniform64_distribution(uint64_t *out);
void uniform64_distribution_vector(uint64_t *out, size_t len);




//uniform distrib between 0 and 1
void random_double(long double *out);
void random_double_vector(long double *out, size_t len);

//Noise distribution
void noise(long double *out, long double param);
void noise_vector(long double *out, long double param, size_t len);

//auxiliary inputs to call discrete Gaussian sampler
gaussian_param_t gaussian(long double param);

//Gaussian over Z
void small_gaussian_overZ(uint64_t *out, long double param);
void small_gaussian_overZ_vector(uint64_t *out, long double param, size_t len);
void gaussian_overZ(uint64_t *out, gaussian_param_t param);
void gaussian_overZ_vector(uint64_t *out, gaussian_param_t param, size_t len);

//G inverse
void gaussian_ginv(uint64_t *out, uint64_t in);
void gaussian_ginv_vector(uint64_t *out, uint64_t *in, size_t len);

//G inverse on polynomials
void gaussian_ginv_poly(uint64_t *out, uint64_t *in);
void gaussian_ginv_poly_vector(uint64_t *out, uint64_t *in, size_t len);

//Precomputes stuff
void precompute_random(gaussian_param_t param);
//Clearing precomputation stuff
void clear_gaussian_param(gaussian_param_t p);
void clear_random();
