#pragma once

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
//#include <complex.h>
//#include <fftw3.h>

//#include "types.h"

//#include "misc.h"

typedef struct
{
	uint64_t* coeffs;
} polynomial_t;

//extern double *r2c_in, *c2r_out, *in_poly_2;
//extern fftw_complex *r2c_out, *c2r_in, *out_poly_1, *out_poly_2;
//extern fftw_plan p_r2c, p_c2r;

void poly_init(polynomial_t *p, size_t size);
void poly_init_zero(polynomial_t *p, size_t size);
void poly_clear(polynomial_t p);
void poly_copy(polynomial_t *dest, polynomial_t src, size_t size);
void poly_copy_index_from(uint64_t *out, uint64_t *temp, size_t index, size_t size);
void poly_copy_index_to(uint64_t *out, uint64_t *temp, size_t index, size_t size);


// scale (multiply or divide) all the coefficients by sc
void poly_scalar_multiply(polynomial_t *out, polynomial_t in, uint64_t *sc, size_t size);
void poly_scalar_divide(polynomial_t  *out, polynomial_t in, uint64_t *sc, size_t size);

void poly_negate(uint64_t *p, size_t size);

//void poly_round(polynomial_t *p, size_t size);

// out = out + in
void poly_add_to(uint64_t *out,uint64_t* in, size_t len);

// out = out - in
void poly_substract_from(polynomial_t *out, polynomial_t in, size_t len);

// out = op1 * op2
//void poly_multiply(polynomial_t *out, polynomial_t op1, polynomial_t op2, size_t N);
void poly_multiply_naive(polynomial_t *out,polynomial_t op1, polynomial_t op2, size_t N);
void poly_multiply_naive_aux(uint64_t *result, uint64_t *op1, uint64_t *op2,size_t N);

// out = out + op1*op2
void poly_add_mult(uint64_t *out,uint64_t* op1,uint64_t* op2,size_t k, size_t len);

// printing function
void print_polynomial(polynomial_t p, size_t len);

//void poly_mul_fft(polynomial_t *out, polynomial_t p1, polynomial_t p2, size_t size);

//void poly_reduce_mod_one(polynomial_t *p, size_t size);

//void poly_decompose(polynomial_t *out, polynomial_t in, size_t size, size_t ell);
