//#pragma once

#include "random.h"

#define _target 3
#define _persample 16

// 20 and 5000 for faster pre-computation
#define _init_fills 20
#define _pool_size 5000
#define _sample_per_fill (_Bg * 250)

//Windows definition for generating random numbers
#ifdef _WIN32
#define RtlGenRandom SystemFunction036 // CSPRNG of windows
BOOLEAN NTAPI RtlGenRandom(PVOID RandomBuffer, ULONG RandomBufferLength); // Used to pull secure random bits on windows
#endif

//Definition of constants
#define TWO_64_MINUS_1 (*(uint64_t*) "\xff\xff\xff\xff\xff\xff\xff\xff")
//long double TWO_64_MINUS_1L = powl(2.0L, 64.0L)-1;

#ifndef M_PI
#define M_PI acosl(-1.0)
#endif

//Structures and global variables for precomputations
struct gaussian_CDT
{
	long double param;
	long double* table;	
	size_t size;
	uint64_t* ss;
	uint64_t* z;
	size_t target;
	uint64_t adjust;
};
struct gaussian_CDT CDT; // pointer table is initialized to NULL

//Precomputed pools of gaussians over Z sorted by residue modulo Bg
struct gaussian_pools
{
	long double param;
	uint64_t mask;
	size_t maxSize;
	size_t* sizes;
	uint64_t** pools;
	//size_t count;
};
struct gaussian_pools Pools; // pointer pools is initialized to NULL
// Hardcoded maxSize in the init function

//uniform distrib over {0,1}
void random_binary(uint64_t *out)
{
	uniform64_distribution(out);
	*out = *out >> 63;
}
void random_binary_vector(uint64_t *out, size_t len)
{
	uniform64_distribution_vector(out, len);
	size_t i;
	for (i=0; i<len; ++i)
		out[i] = out[i] >> 63;
}

//uniform distrib between 0 and 2^64-1
void uniform64_distribution(uint64_t *out)
{
	#ifdef _WIN32
	RtlGenRandom((PVOID) out,(ULONG) sizeof(uint64_t));
    #else// __unix__ __macos__
    int fd=open("/dev/urandom", O_RDONLY);
	
	read(fd,out,sizeof(uint64_t));
	
	close(fd);
    #endif
}

void uniform64_distribution_vector(uint64_t *out, size_t len)
{
	#ifdef _WIN32
	RtlGenRandom((PVOID) out,(ULONG) len*sizeof(uint64_t));
	#else// __unix__ __macos__
	int fd=open("/dev/urandom", O_RDONLY);
	read(fd,out,len*sizeof(uint64_t));
	close(fd);
    #endif
}


//uniform distrib between 0 and 1
void random_double(long double *out)
{
	uint64_t tmp;
	uniform64_distribution(&tmp);
	*out = tmp / (long double) TWO_64_MINUS_1;
}
void random_double_vector(long double *out, size_t len)
{
	uint64_t *tmp = (uint64_t*) malloc(len * sizeof(uint64_t));
	uniform64_distribution_vector(tmp, len);
	size_t i;
	for (i=0; i<len; ++i)
		out[i] = tmp[i] / (long double) TWO_64_MINUS_1;
	free(tmp);
}

//Noise distribution
void noise(long double *out, long double param)
{
	uint64_t coin;
	random_binary(&coin);
	long double tmp;
	random_double(&tmp);
	if (coin)
		tmp = param * sqrtl(-2. * logl(tmp));
	else
		tmp = - param *sqrtl(-2. * logl(tmp));
	*out = tmp;
}


void noise_vector(long double *out, long double param, size_t len)
{
	long double *tmp = (long double*) malloc(len * sizeof(long double));
	random_double_vector(tmp, len);
	size_t i;
	for (i=0; i<len-1; i+=2) {
	  tmp[i] = sqrtl(-2. * logl(tmp[i]));
	  tmp[i+1] = 2. * M_PI * tmp[i+1];
	  out[i] = (param * tmp[i] * cosl(tmp[i+1]));
	  out[i+1] = param * tmp[i] * sinl(tmp[i+1]);
	}
	free(tmp);
	if (len % 2)
		noise(out+len-1, param);
}

//Gaussian over Z
void precompute_CDT(long double param)
{
	CDT.param = param;
	CDT.size = 2048;
	CDT.table = (long double*) malloc (CDT.size * sizeof(long double));
	size_t i;
	long double norm = 0.0L;

	for (i=0; i<CDT.size; ++i)
	{
		norm += expl(- M_PI * powl(((long double) ((int64_t) (i-CDT.size/2))/64),2.0L));
		CDT.table[i] = norm;
	}

	for (i=0; i<CDT.size; ++i)
	{
		CDT.table[i] /= norm;
//			printf("%Id : %f\n", i-CDT.size/2, (double) CDT.table[i]);   // Uncomment to print cumulative distribution table
	}

	CDT.target = _target;

	CDT.ss = (uint64_t*) malloc (CDT.target * sizeof(uint64_t));
	CDT.z = (uint64_t*) malloc (CDT.target * sizeof(uint64_t));
	uint64_t ss=4096,z;
	for (i=0; i<CDT.target; ++i)
	{
		z = (sqrtl(ss)/16.0L);
		CDT.z[i]=z;
		ss = (2*z*(z-1)+1)*ss;
		CDT.ss[i]=ss;
//		printf("z=%llu,ss=%llu\n",z,ss);
	}
	CDT.adjust = 1+(1.0L+sqrtl(2*param*param/CDT.ss[CDT.target-1]-1))/2.0L;
}

void clear_CDT()
{
	if (CDT.table != NULL)
	{
		free(CDT.table);
		CDT.table = NULL;
		free(CDT.ss);
		free(CDT.z);
	}
}

void small_gaussian_overZ(uint64_t *out, long double param)
{
	if (CDT.table == NULL) //Cumulative distribution table has not yet been initialized, so we do it
	{
		precompute_CDT(param);
	}

	long double coin;
	random_double(&coin);

	size_t a = 0;
  	size_t b = CDT.size-1;
  	size_t c;
  	while (b-a > 1)
  	{
  		c = (a+b)/2;
  		if (coin > CDT.table[c])
  			a = c;
  		else
  			b = c;
  	}
  	*out = b-CDT.size/2;
}

void small_gaussian_overZ_vector(uint64_t *out, long double param, size_t len)
{
	if (CDT.table == NULL) //Cumulative distribution table has not yet been initialized, so we do it
	{
		precompute_CDT(param);
	}

	long double *coin = (long double*) malloc (len * sizeof(long double));
	random_double_vector(coin, len);

	size_t i;
	for (i=0; i<len; ++i)
	{
		size_t a = 0;
  		size_t b = CDT.size-1;
  		size_t c;
  		while (b-a > 1)
  		{
  			c = (a+b)/2;
  			if (coin[i] > CDT.table[c])
  				a = c;
  			else
  				b = c;
  		}
  		out[i] = b-CDT.size/2;
	}
	free(coin);
}

uint64_t combine(uint64_t* in,  size_t size, size_t level)
{
	if (!level)
		return *in;
	uint64_t l = level-1;
	uint64_t z = CDT.z[l];
	uint64_t newsize = size/2;
	return z*combine(in, newsize, l)+(z-1)*combine(in+newsize, newsize, l);
}

void gaussian_overZ(uint64_t *out, long double param)
{
	uint64_t* sample = (uint64_t *) malloc(_persample*sizeof(uint64_t));
	small_gaussian_overZ_vector(sample, param, _persample);
	*out = CDT.adjust * combine(sample,8,3) + (CDT.adjust-1) * combine(sample+8,8,3);
	free(sample);
}

void gaussian_overZ_vector(uint64_t *out, long double param, size_t len)
{
	uint64_t* sample = (uint64_t *) malloc(_persample*len*sizeof(uint64_t));
	small_gaussian_overZ_vector(sample, param, _persample*len);
	for (size_t i=0;i<len;++i)
		out[i] = CDT.adjust * combine(sample+16*i,8,3) + (CDT.adjust-1) * combine(sample+16*i+8,8,3);
	free(sample);
}

//G inverse
/* G SPECS */
/* (d+1)ell ROWS, (d+1)N COLUMNS, Each adjacent N numbers are polynomial coefficients, Powers of Bg decreasing top to bottom from Bg^(ell-1) to 1 */
void init_pools(long double param)
{
	printf("Initializing pools\n");
	Pools.maxSize = _pool_size;
	Pools.param = param;
	Pools.mask = (1 << _logBg) - 1;
	Pools.sizes = (size_t *) malloc(_Bg * sizeof(size_t));
	Pools.pools = (uint64_t **) malloc(_Bg * sizeof(uint64_t *));
	//Pools.count = 0;
	size_t i;
	for (i=0; i<_Bg; ++i)
	{
		Pools.sizes[i] = 0;
		Pools.pools[i] = (uint64_t *) malloc(Pools.maxSize * sizeof(uint64_t));
	}
	printf("Pools initialized\n");
}

void clear_pools()
{
	if (Pools.pools != NULL)
	{
		size_t i;
		for (i = 0; i < _Bg; ++i)
		{
			free(Pools.pools[i]);
		}	
		free(Pools.pools);
		free(Pools.sizes);
		Pools.pools = NULL;
		Pools.sizes = NULL;	
	}
}

void fill_pool(long double param) // Fills with _sample_per_fill samples the pools
{
	//printf("Filling pools\n");
	uint64_t* sample = (uint64_t *) malloc(_sample_per_fill*sizeof(uint64_t));
	gaussian_overZ_vector(sample, param, _sample_per_fill);
	uint64_t v;
	size_t i;
	for (i=0; i<_sample_per_fill; ++i)
	{
		v = sample[i] % _Bg;
		if (Pools.sizes[v] < Pools.maxSize - 1)
		{
			Pools.pools[v][++Pools.sizes[v]] = sample[i];
			//++Pools.sizes[t];
		}
	}
	free(sample);
/*	printf("Pools filled :\n");
	for (t=0; t<Bg; ++t)
	{
		printf("%Iu : %Iu elements\n",t,Pools.sizes[t]);
		for (i=1; i < 6; ++i)
		{
			if (i<Pools.sizes[t])
				printf("\t%Id\n",Pools.pools[t][Pools.sizes[t]-i]);
		}
	}*/
	//Pools.count++;
	//printf("%I64u\n", Pools.count);
}

void pull_pool(uint64_t *out, uint64_t in, long double param)
{
	//printf("%I64u\n", in);
	while (!Pools.sizes[in])
	{
		fill_pool(param);
	}
	*out = Pools.pools[in][--Pools.sizes[in]]; // This line for real precomputation
	//*out = Pools.pools[in][Pools.sizes[in]-1]; // This line for dummy precomputation
	//--Pools.sizes[in];
	//Pools.count++;
	//printf("%I64u\n", Pools.count);
}

void gaussian_ginv(uint64_t *out, uint64_t in, long double param)
{
	if (Pools.pools == NULL) //Sample pools have not yet been initialized, so we do it
	{
		init_pools(param);
	}

	
	//printf("v : %Iu\n", v);
	size_t i;
	/*for (i=ell-1; i+1 > 0; --i)
	{
		//printf("v mod Bg : %Id\n", v%Bg);
		pull_pool(out+i, v%Bg, param, Bg);
		//pull_pool(out+i, (in << (i * Pools.logBg)) >> (64 - Pools.logBg), param, Bg);
		//v = (v - out[i])/((int64_t) Bg);
		v = (v - out[i])/Bg;
		//printf("v : %Id\n", v);
	}*/
	uint64_t v = in >> (64 - _logBg * _ell);
	for (i = _ell-1; i+1 > 0; --i)
	{
		pull_pool(out+i, v & Pools.mask, param);
		v = (v - out[i]) >> _logBg;
	}
}

void gaussian_ginv_vector(uint64_t *out, uint64_t *in, long double param, size_t len)
{
	size_t i;
	for(i=0; i<len; ++i)
		gaussian_ginv(out+(i*_ell), in[i], param);
}

//N first elements of output is LSBs
void gaussian_ginv_poly(uint64_t *out, uint64_t *in, long double param)
{
	if (Pools.pools == NULL) //Sample pools have not yet been initialized, so we do it
	{
		init_pools(param);
	}
	size_t coeff;
	size_t i;
	for (coeff = 0; coeff < _N; ++coeff)
	{
		int64_t v = in[coeff] >> (64 - _logBg * _ell);
		for (i=_ell-1; i+1 > 0; --i)
		{
			//pull_pool(out+i*N+coeff, v%Bg, param, Bg);
			pull_pool(out+i*_N+coeff, v & Pools.mask, param);
			v = (v - out[i*_N+coeff]) >> _logBg;
			//v = (v - out[i*N+coeff])/Bg;
		}
	}
}

void gaussian_ginv_poly_vector(uint64_t *out, uint64_t *in, long double param, size_t len)
{
	for(size_t i = 0; i < len; ++i)
		gaussian_ginv_poly(out+(i*_ell*_N), in + i*_N, param);
}


//Initialize and fills pools for Ginv
void precompute_random(long double param)
{
	init_pools(param);
	for (size_t i = 0; i < _init_fills; ++i)
		fill_pool(param);
}

// Clearing the precomputation stuff
void clear_random()
{
	clear_CDT();
	clear_pools();
}