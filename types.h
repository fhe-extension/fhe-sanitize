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
typedef tlwe_sample *pkg;//m encryptions of 0
//typedef tlwe_sample epk;//Linear combination of m encryptions of 0 // Same as tlwe_sample
