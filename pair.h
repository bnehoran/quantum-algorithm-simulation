#ifndef PAIR_H
#define PAIR_H

#include <complex.h>
#include "bignum.h"

/* A Pair is an (index,value) pair used in the sparse 
   vector representation */

struct Pair {
    long index;
    bignum index_bn;
    double complex val;
};

#endif