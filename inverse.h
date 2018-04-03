#ifndef INVERSE_H
#define INVERSE_H

#include "operator.h"
#include "state.h"
#include "vector.h"
#include "pair.h"
#include <complex.h>
#include <stdio.h>

/* Returns an initial approximation for 1/w that is 2^-p for p such that
   2^p > w >= 2^(p-1).
   w is represented as a fixed precision number stored in the wr Register
   of State s. xr is the Register of s where the result (2^-p) is stored and 
   anc is the index of some ancilla qubit (the qubits represented by xr and 
   anc must all be initialized to 0) */
State Inverse_init(State s, Register wr, Register xr, int anc);

/* Takes an positive integer w, encodes it as a quantum mechanical state,
   and performs Newton's iterations to generate successively better 
   approximations for 1/w.
   w is represented as a fixed precision number with m integer bits and
   n bits total.
   b is the desired number of fractional bits in the result. b >= n-m. */
State Inverse_inverse(long w, int n, int m, int b);


void Inverse_run_tests(void);


#endif