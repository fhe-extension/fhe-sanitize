//#pragma once

#include "random.h"
#define _target_max 4
#define _tail_cut 6

#define _threads_number 2
#ifdef _threads_number
#include <pthread.h>
#include <semaphore.h>
#endif

#define _buffer_size 25600000
#define _sample_size 25600000

#define _init_fills 2
#define _pool_size 5000
#define _sample_per_fill (_Bg * 2500)

//Windows definition for generating random numbers
#ifdef _WIN32
#define RtlGenRandom SystemFunction036 // CSPRNG of windows
BOOLEAN NTAPI RtlGenRandom(PVOID RandomBuffer, ULONG RandomBufferLength); // Used to pull secure random bits on windows
#endif

//Definition of constants
#define TWO_64_MINUS_1 (*(uint64_t*) "\xff\xff\xff\xff\xff\xff\xff\xff")
//double TWO_64_MINUS_1L = powl(2.0L, 64.0L)-1;

#ifndef M_PI
#define M_PI acosl(-1.0)
#endif

#ifdef _threads_number
pthread_t threads[_threads_number - 1];
#endif

//Structures and global variables for precomputations
struct gaussian_CDT
{
  double param;
  double* table;
  int size;
};
struct gaussian_CDT CDT; // pointer table is initialized to NULL

//Precomputed pools of gaussians over Z sorted by residue modulo Bg
struct gaussian_pools
{
  gaussian_param_t param;
  uint64_t mask;
  int maxSize;
  int* sizes;
  uint64_t** pools;
#ifdef _threads_number
  pthread_mutex_t *mutex;
  //sem_t *fill;
  sem_t *pull;
#endif
};
struct gaussian_pools Pools; // pointer pools is initialized to NULL
// Hardcoded maxSize in the init function

double *coin = NULL;
uint64_t *sample = NULL;
uint64_t *buffer = NULL;

//uniform distrib over {0,1}
void random_binary_vector(uint64_t *out, int len)
{
  uniform64_distribution_vector(out, len);
  for (int i = 0; i < len; ++i)
    out[i] = out[i] >> 63;
}
void random_binary(uint64_t *out)
{
  random_binary_vector(out, 1);
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

void uniform64_distribution_vector(uint64_t *out, int len)
{
#ifdef _WIN32
  RtlGenRandom((PVOID) out,(ULONG) len*sizeof(uint64_t));
#else// __unix__ __macos__
  //int fd=open("/dev/urandom", O_RDONLY);
  //read(fd,out,len*sizeof(uint64_t));
  //close(fd);
  //getrandom(out, len * sizeof(uint64_t), GRND_NONBLOCK);
  RAND_bytes((void *)out, len * sizeof(uint64_t));
#endif
}


//uniform distrib between 0 and 1
void random_double(double *out)
{
  uint64_t tmp;
  uniform64_distribution(&tmp);
  *out = tmp / (double) TWO_64_MINUS_1;
}
void random_double_vector(double *out, int len)
{
  //printf("%ld random_double_vector\n", len);
  if (!buffer)
    buffer = malloc(_buffer_size * sizeof(uint64_t));
  uniform64_distribution_vector(buffer, len);
  for (int i = 0; i < len; ++i)
    out[i] = buffer[i] / (double) TWO_64_MINUS_1;
}

//Noise distribution
void noise(double *out, double param)
{
  uint64_t coin;
  random_binary(&coin);
  double tmp;
  random_double(&tmp);
  if (coin)
    tmp = param * sqrtl(-2. * logl(tmp));
  else
    tmp = - param *sqrtl(-2. * logl(tmp));
  *out = tmp;
}


void noise_vector(double *out, double param, int len)
{
  double *tmp = (double*) malloc(len * sizeof(double));
  random_double_vector(tmp, len);
  for (int i = 0; i < len - 1; i += 2) {
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
void check_CDT(double *table, int size){
  double expectation = 0.0;
  for (int i = 0; i < size-1; i = i + 1){
    //printf("Pr(%d) = %lf\n", i - size/2, table[i+1] - table[i]);
    expectation += (i - size/2) * (table[i+1] - table[i]);
  }
  printf("AVG : %lf\n", expectation);
  double variance = 0.0;
  for (int i = 0; i < size-1; i = i + 1){
    variance += (i - size/2 - expectation) * (i - size/2 - expectation) * (table[i+1] - table[i]);
  }
  printf("VAR : %lf\n", variance);
}

double *compute_CDT(double param)
{
  //printf("compute_CDT\n");
  int size = 2*_tail_cut*param;
  double *table = (double*) malloc (size * sizeof(double));
  double norm = 0.0L;

  for (int i = 0; i < size; ++i)
    {
      table[i] = norm;
      norm += expl(- M_PI * powl(((double) ((int64_t) (i-size/2))/(param)),2.0L));
    }

  for (int i = 0; i < size; ++i)
    {
      table[i] /= norm;
      //printf("%Id : %f\n", i-size/2, (double) table[i]);   // Uncomment to print cumulative distribution table
    }
  /*printf("Param : %lf\n",param);
  printf("Expected variance : %lf\n",param * param / (2*acos(-1.0L)));
  check_CDT(table, size);*/
  return table;
}

void precompute_CDT(double param) {
  CDT.param = param;
  CDT.size = 2*_tail_cut*param;
  CDT.table = compute_CDT(param); //(double*) malloc (CDT.size * sizeof(double));
}

gaussian_param_t gaussian(double param) {
  if (param < 250)
    return (gaussian_param_t) {
      .param = param,
      .target = 0,
      .z = NULL,
      .adjust = 1,
    };

  double ssi, next;
  uint64_t zi; int t;

  //printf("Param : %lf - squared : %lf\n", param, param * param);

  double best_s0 = 34;
  double best_tildes = 1.0 / 0.0;
  for(double s0 = 34; s0 <= 512; s0 += 0.01){
    t = 0;
    ssi = s0*s0;
    zi = (sqrtl(ssi)/16.0L);
    next = (2*zi*(zi-1)+1)*ssi;
    //printf("Test s0 = %lf - ",s0);
    while (param*param > next) {
      t++;
      ssi = next;
      zi = floorl((sqrtl(ssi)/16.0L));
      next = (zi > 2 ? 2*zi*(zi - 1) + 1 : zi * zi + 1)*ssi;
    }
    uint64_t adjustz = ceill(1+(1.0L+sqrtl(2*param*param/ssi-1))/2.0L);
    double tildes = (2 * adjustz * (adjustz - 1) + 1) * ssi;
    if (tildes < best_tildes){
      best_tildes = tildes;
      best_s0 = s0;
    }
    //printf("tildes : %lf\n", tildes);
    //printf("OK range : %lf -> %lf, next is %lf\n", ssi, 5 * ssi, next);
  }
  gaussian_param_t p;
  p.param = best_s0;
  p.target = t;
  //printf("s0=%lf\n", p.param);
  p.z = malloc(t*sizeof(uint64_t));
  ssi = best_s0*best_s0;
  //printf("ss0=%lf\n", ssi);
  for (int i = 0; i < t; ++i) {
    zi = floorl((sqrtl(ssi)/16.0L));
    p.z[i] = zi;
    ssi = (zi > 2 ? 2*zi*(zi - 1) + 1 : zi * zi + 1)*ssi;
    //printf("z%d=%lu,ss%d=%lf\n",1+i,zi,1+i,ssi);
  }
  p.adjust = ceill(1+(1.0L+sqrtl(2*param*param/ssi-1))/2.0L);
  //printf("adjust z = %ld\n", p.adjust);
  return p;
}

void clear_CDT()
{
  if (CDT.table != NULL)
    {
      free(CDT.table);
      CDT.table = NULL;
    }
}

void small_gaussian_overZ(uint64_t *out, double param)
{
  small_gaussian_overZ_vector(out, param, 1);
}

void small_gaussian_overZ_vector(uint64_t *out, double param, int len)
{
  double *table;
  int size;
  if (CDT.param == param) {//We want to sample from the precomputed gaussian parameter
    table = CDT.table;
    size = CDT.size;
    //printf("%ld\n", size);
  }
  else {
    table = compute_CDT(param);
    size = 2*_tail_cut*param;
  }

  if (!coin)
    coin = (double*) malloc (_sample_size * sizeof(double));
  random_double_vector(coin, len);

  for (int i = 0; i < len; ++i) {
    uint64_t start = 0;
    uint64_t len = size;
    while (len > 2)
      {
        if (coin[i] >= table[start+len/2]){
          start += len / 2;
          len = len - len / 2;
        }
        else
          len = len / 2;
      }
    out[i] = start-size/2;
  }

  if (CDT.param != param)
    free(table);
}

uint64_t combine(uint64_t* in, int size, int level, uint64_t *z)
{
  if (!level)
    return *in;
  uint64_t l = level-1;
  uint64_t newsize = size/2;
  return z[l]*combine(in, newsize, l, z)+(z[l]-1)*combine(in+newsize, newsize, l, z);
}

void gaussian_overZ(uint64_t *out, gaussian_param_t param)
{
  int n_samples = 1 << param.target;
  uint64_t* sample = (uint64_t *) malloc(2*n_samples*sizeof(uint64_t));
  small_gaussian_overZ_vector(sample, param.param, 2*n_samples);
  *out = param.adjust * combine(sample,n_samples,param.target,param.z) + (param.adjust-1) * combine(sample+n_samples,n_samples,param.target,param.z);
  free(sample);
}

void gaussian_overZ_vector(uint64_t *out, gaussian_param_t param, int len)
{
  //Pools.count++;
  //printf("%I64u\n", Pools.count);
  int n_samples = 1 << param.target;
  if (!sample)
    sample = (uint64_t *) malloc(_sample_size*sizeof(uint64_t));
  small_gaussian_overZ_vector(sample, param.param, 2*n_samples*len);

  for (int l = 0; l < param.target; ++l){
    for (int i = 0; i < n_samples * len; ++i)
      sample[i] = param.z[l] * sample[i] + (param.z[l] > 2 ? param.z[l] - 1 : 1) * sample[i + n_samples * len];
    n_samples >>= 1;
  }
  for (int i = 0; i < len; ++i)
    out[i] = param.adjust * sample[i] + (param.adjust - 1) * sample[i + len];
}

//G inverse
/* G SPECS */
/* (d+1)ell ROWS, (d+1)N COLUMNS, Each adjacent N numbers are polynomial coefficients, Powers of Bg decreasing top to bottom from Bg^(ell-1) to 1 */
void init_pools(gaussian_param_t param)
{
  printf("Initializing pools\n");
  precompute_CDT(param.param);
  Pools.maxSize = _pool_size;
  Pools.param = param;
  Pools.mask = (1 << _logBg) - 1;
  Pools.sizes = (int *) malloc(_Bg * sizeof(int));
  Pools.pools = (uint64_t **) malloc(_Bg * sizeof(uint64_t *));
#ifdef _threads_number
  Pools.mutex = malloc(_Bg * sizeof(pthread_mutex_t));
  //Pools.fill = malloc(_Bg * sizeof(sem_t));
  Pools.pull = malloc(_Bg * sizeof(sem_t));
#endif
  for (int i = 0; i < _Bg; ++i)
    {
      Pools.sizes[i] = 0;
      Pools.pools[i] = (uint64_t *) malloc(Pools.maxSize * sizeof(uint64_t));
#ifdef _threads_number
      pthread_mutex_init(Pools.mutex + i, NULL);
      //sem_init(Pools.fill + i, 0, Pools.maxSize);
      sem_init(Pools.pull + i, 0, 0);
#endif
    }
  printf("Pools initialized\n");
}

void clear_gaussian_param(gaussian_param_t p) {
  if (p.z != NULL)
    free(p.z);
}

void clear_pools()
{
  if (Pools.pools != NULL)
    {
      for (int i = 0; i < _Bg; ++i)
        {
          free(Pools.pools[i]);
        }
      free(Pools.pools);
      free(Pools.sizes);
      Pools.pools = NULL;
      Pools.sizes = NULL;
      clear_gaussian_param(Pools.param);
    }
}

void *fill_pool(void *args) // Fills with _sample_per_fill samples the pools
{
  assert(!args);
  //printf("Filling pools\n");
  uint64_t* sample = (uint64_t *) malloc(_sample_per_fill*sizeof(uint64_t));
#ifdef _threads_number
  for (;;){
#endif
  gaussian_overZ_vector(sample, Pools.param, _sample_per_fill);
  uint64_t v;
  for (int i = 0; i < _sample_per_fill; ++i)
    {
      v = sample[i] % _Bg;
      //assert(v >= 0 && v < _Bg);
      //printf("SAMPLED : %lu\n", v);
#ifdef _threads_number
      //sem_wait(Pools.fill + v);
      pthread_mutex_lock(Pools.mutex + v);
      if (Pools.sizes[v] < Pools.maxSize){
        Pools.pools[v][Pools.sizes[v]++] = sample[i];
        sem_post(Pools.pull + v);
      }
      pthread_mutex_unlock(Pools.mutex + v);
#else
      if (Pools.sizes[v] < Pools.maxSize)
        {
          Pools.pools[v][Pools.sizes[v]++] = sample[i];
        }
#endif
    }
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
#ifdef _threads_number
  }
#endif
  free(sample);
  return NULL;
}

void pull_pool(uint64_t *out, uint64_t in)
{
  //assert(0 <= in && in < _Bg);
#ifdef _threads_number
  sem_wait(Pools.pull + in);
  pthread_mutex_lock(Pools.mutex + in);
  *out = Pools.pools[in][--Pools.sizes[in]];
  pthread_mutex_unlock(Pools.mutex + in);
  //sem_post(Pools.fill + in);
#else
  while (!Pools.sizes[in])
    {
      fill_pool(NULL);
    }
  *out = Pools.pools[in][--Pools.sizes[in]]; // This line for real precomputation
  //*out = Pools.pools[in][Pools.sizes[in]-1]; // This line for dummy precomputation
#endif
}

void gaussian_ginv(uint64_t *out, uint64_t in)
{
  //	if (Pools.pools == NULL) //Sample pools have not yet been initialized, so we do it
  //	{
  //		init_pools(param);
  //	}
  uint64_t v = in >> (64 - _logBg * _ell);
  for (int i = _ell-1; i+1 > 0; --i)
    {
      pull_pool(out+i, v & Pools.mask);
      v = (v - out[i]) >> _logBg;
    }
}

void gaussian_ginv_vector(uint64_t *out, uint64_t *in, int len)
{
  for(int i = 0; i < len; ++i)
    gaussian_ginv(out+(i*_ell), in[i]);
}

//N first elements of output is LSBs
void gaussian_ginv_poly(uint64_t *out, uint64_t *in)
{
  //	if (Pools.pools == NULL) //Sample pools have not yet been initialized, so we do it
  //	{
  //		init_pools(param);
  //	}
  for (int coeff = 0; coeff < _N; ++coeff)
    {
      int64_t v = in[coeff] >> (64 - _logBg * _ell);
      for (int i = _ell-1; i + 1 > 0; --i)
        {
          //printf("PULL : %lu\n", v & Pools.mask);
          pull_pool(out+i*_N+coeff, v & Pools.mask);
          //printf("GOT\n");
          v = (v - out[i*_N+coeff]) >> _logBg;
        }
    }
}

void gaussian_ginv_poly_vector(uint64_t *out, uint64_t *in, int len)
{
  for(int i = 0; i < len; ++i)
    gaussian_ginv_poly(out+(i*_ell*_N), in + i*_N);
}


//Initialize and fills pools for Ginv
void precompute_random(gaussian_param_t param)
{
  init_pools(param);
#ifdef _threads_number
  for (int i = 0; i < _threads_number - 1; ++i){
    pthread_create(threads + i, NULL, fill_pool, NULL);
  }
#else
  for (int i = 0; i < _init_fills; ++i)
    fill_pool(NULL);
#endif
}

// Clearing the precomputation stuff
void clear_random()
{
  clear_CDT();
  clear_pools();
  if (coin) free(coin);
  if (sample) free(sample);
  if (buffer) free(buffer);
}
