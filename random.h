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
  double param;
  size_t target;
  uint64_t* z;
  uint64_t adjust;
};
typedef struct gaussian_param gaussian_param_t;

//uniform distrib over {0,1}
void random_binary(uint64_t *out);
void random_binary_vector(uint64_t *out, size_t len);

//uniform distrib between 0 and 2^64-1
void uniform64_distribution(uint64_t *out);
void uniform64_distribution_vector(uint64_t *out, size_t len);

//uniform distrib between 0 and 1
void random_double(double *out);
void random_double_vector(double *out, size_t len);

//Noise distribution
void noise(double *out, double param);
void noise_vector(double *out, double param, size_t len);

//auxiliary inputs to call discrete Gaussian sampler
gaussian_param_t gaussian(double param);

//Gaussian over Z
void small_gaussian_overZ(uint64_t *out, double param);
void small_gaussian_overZ_vector(uint64_t *out, double param, size_t len);
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
