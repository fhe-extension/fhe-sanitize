//#pragma once

#include "polynomial.h"

#define B_g 1024
#define ELL 3
#define PREC 1073741824
#define TWOPREC 2147483648
#define LOGB 10


// declare global variables for FFT
//double *r2c_in, *c2r_out, *in_poly_1, *in_poly_2;
//fftw_complex *r2c_out, *c2r_in, *out_poly_1, *out_poly_2;
//fftw_plan p_r2c, p_c2r;


void poly_init(polynomial_t *p, size_t size)
{
	p->coeffs = (uint64_t*) malloc(size * sizeof(uint64_t));
}

void poly_init_zero(polynomial_t *p, size_t size)
{
	p->coeffs = (uint64_t*) malloc(size * sizeof(uint64_t));
	for(size_t i = 0; i < size; ++i) {
		p->coeffs[i] = 0;
	}
}

void poly_clear(polynomial_t p)
{
	free(p.coeffs);
}

void poly_copy(polynomial_t *out, polynomial_t temp, size_t size)
{
	for(size_t i = 0; i < size; ++i) {
		out->coeffs[i] = temp.coeffs[i];
	}
}


void poly_copy_index_from(uint64_t *out, uint64_t *temp, size_t index, size_t size)
{
	for(size_t i = 0; i < size; ++i) {
		out[index+i] = temp[i];
	}
}

void poly_copy_index_to(uint64_t *out, uint64_t *temp, size_t index, size_t size)
{
	for(size_t i = 0; i < size; ++i) {
		out[i] = temp[index+i];
	}
}



void poly_scalar_multiply(polynomial_t *out, polynomial_t in, uint64_t* sc, size_t size)
{
	for(size_t i = 0; i < size; ++i) {
		out->coeffs[i] = in.coeffs[i]* (*sc);
	}
}

void poly_scalar_divide(polynomial_t *out, polynomial_t in, uint64_t* sc, size_t size)
{
	for(size_t i = 0; i < size; ++i) {
		out->coeffs[i] = in.coeffs[i] / (*sc);
	}
}

/*void poly_negate(polynomial_t *p, size_t size)
{
	for (size_t i = 0; i < size; ++i) {
		p->coeffs[i] = -(p->coeffs[i]);
	}
}
*/


void poly_add_to(uint64_t *out, uint64_t* in, size_t len)
{
	for(size_t i = 0; i < len; ++i) {
		out[i] += in[i];
	}
}

void poly_sub_from(polynomial_t *out, polynomial_t in,size_t len)
{

	for(size_t i = 0; i < len; ++i) {
		out->coeffs[i] -= in.coeffs[i];
	}
}

//void poly_mul(polynomial_t *out, polynomial_t op1, polynomial_t op2, size_t N)
//{
//	poly_mul_fft(out, op1, op2, N);
//}

void poly_multiply_naive(polynomial_t *out, polynomial_t op1, polynomial_t op2,size_t N)
{
	
	int degree, sign;
	sign = 1;
	for(size_t i = 0; i < N; ++i)
		out->coeffs[i] = 0;
	for(size_t i = 0; i < N; ++i) {
		for(size_t j = 0; j < N; ++j) {
			if(i + j >= N) {
				degree = i + j - N;
				sign = -1;
			}
			else {
				degree = i + j;
				sign = 1;
			}
			out->coeffs[degree] += sign * (op1.coeffs[i] * op2.coeffs[j]);
		}
	}
	
}

void poly_multiply_naive_aux(uint64_t *result, uint64_t *op1,uint64_t *op2,size_t N) {
	uint64_t ri;
    for (size_t i=0; i<N; i++) {
		ri=0;
			for (size_t j=0; j<=i; j++) {
		    	ri += op1[j]*op2[i-j];
			}
			for (size_t j=i+1; j<N; j++) {
		    	ri -= op1[j]*op2[N+i-j];
			}
		result[i]=ri;
    }
}



/*void poly_add_mult(uint64_t* out, uint64_t* op1, uint64_t* op2, size_t k,size_t len)
{
	polynomial_t tmp;
	poly_init(&tmp, len);
	poly_multiply_naive(tmp.coeffs, op1, op2,k,len);
	poly_add_to(out,tmp.coeffs, len);
	poly_clear(tmp);
}*/

void print_polynomial(polynomial_t p,size_t len)
{
	for(size_t i = 0; i < len; ++i) {
		printf("%llu  ", p.coeffs[i]);
	}
	printf("\n");
}











