#pragma once

//#include <stdlib.h>

//#include "random.h"

/* LWE types */
typedef uint64_t *lwe_sample;
typedef uint64_t *lwe_sk;

/* TLWE types */
typedef uint64_t *tlwe_sample;
typedef uint64_t *tlwe_sk;

/* TGSW types */
typedef uint64_t *tgsw_sample;

/* Bootstrapping keys types */
typedef lwe_sample *ksk;
typedef tgsw_sample *bsk;

/* Sanitization keys */
typedef tlwe_sample *pks; // public keys for statistical rerandomization
typedef uint64_t *pkc; // public keys for computational rerandomization
