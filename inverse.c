/* 
 * inverse.c
 * --------
 * A sample client to the library demonstrating the implementation of the
 * inversion algorithm from "Quantum Algorithms and Circuits for Scientific 
 * Computing" by Mihir K. Bhaskar, Stuart Hadfield, Anargyros Papageorgiou,
 * and Iasonas Petras (listed as Algorithm 0), which can be found at 
 * https://arxiv.org/pdf/1511.08253.pdf
 * 
 * 
 * Author: Barak Nehoran
 * Advisor: Iasonas Petras
 */


#include "inverse.h"
#include "operator.h"
#include "state.h"
#include "vector.h"
#include <stdlib.h>
#include <assert.h>
#include <complex.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <string.h>


#define VERBOSE false

int wasted_qubits = 0;

int num_qubits(int n, int m, int b) {
   return m + b + 1 + b;
}

State apply(Operator op, State s) {
   if (VERBOSE)
      return Operator_print_apply(op, s);
   else 
      return Operator_apply(op, s);
}

/* Returns an initial approximation for 1/w that is 2^-p for p such that
   2^p > w >= 2^(p-1).
   w is represented as a fixed precision number stored in the wr Register
   of State s. xr is the Register of s where the result (2^-p) is stored and 
   anc is the index of some ancilla qubit (the qubits represented by xr and 
   anc must all be initialized to 0) */
State Inverse_init(State s, Register wr, Register xr, int anc) {
   assert (Register_isValid(wr));
   assert (Register_isValid(xr));
   assert (wr->high <= State_num_qubits(s) && xr->low <= State_num_qubits(s));
   assert (Register_have_no_overlap(wr,xr)); // check for overlap
   assert (anc >= wr->high || anc < wr->low);
   assert (anc >= xr->high || anc < xr->low);
   int num_q = State_num_qubits(s);
   int m = wr->high - wr->zero;
   int x_b = xr->low;
   int x_1 = xr->zero - 1;
   int x_m = xr->zero - m;
   int w_m_1 = wr->high - 1;
   int w_b = wr->low;
   int w_0 = wr->zero;
   
   
   Operator op;
   
   int cont_1[m+3];
   int cont_0[m+3];
   int targ[m+3];
   cont_1[0] = -1;
   cont_1[1] = -1;
   cont_1[2] = -1;
   cont_0[0] = -1;
   cont_0[1] = -1;
   cont_0[2] = -1;
   targ[0] = -1;
   targ[1] = -1;
   targ[2] = -1;
   int i = w_m_1;
   int j = x_m;
   for (; j < anc; i--, j++) {
      cont_1[0] = i;
      cont_1[1] = -1;
      cont_0[0] = anc;
      cont_0[1] = -1;
      targ[0] = j;
      targ[1] = -1;
      targ[2] = -1;
      op = Operator_cnot_0 (num_q, cont_1, cont_0, targ);
      s = apply(op, s);
      Operator_free(op);
      
      cont_1[0] = j;
      cont_1[1] = -1;
      cont_0[0] = -1;
      cont_0[1] = -1;
      targ[0] = anc;
      targ[1] = -1;
      targ[2] = -1;
      op = Operator_cnot_0(num_q, cont_1, cont_0, targ);
      s = apply(op, s);
      Operator_free(op);
   }
   i = 0;
   j = x_m;
   for (; j < anc; i++, j++) {
      cont_0[i] = j;
   }
   cont_0[i] = -1;
   cont_1[0] = -1;
   targ[0] = anc;
   targ[1] = -1;
   op = Operator_cnot_0(num_q, cont_1, cont_0, targ);
   s = apply(op, s);
   Operator_free(op);
   
   op = Operator_pauli_X(num_q, anc);
   s = apply(op, s);
   Operator_free(op);
   
   return s;
}


/* Takes an positive integer w, encodes it as a quantum mechanical state,
   and performs Newton's iterations to generate successively better 
   approximations for 1/w.
   w is represented as a fixed precision number with m integer bits and
   n bits total. w, must have a value grater than or equal to 1, so that 
   its inverse is fractional.
   b is the desired number of fractional bits in the result. b >= n-m. */
State Inverse_inverse(long w, int n, int m, int b) {
   assert (w >= 0 && n >= 0 && m >= 0 && b >= 0);
   assert (b >= n-m && b >= m);
   assert (w >> (n-m) >= 1L);
   
   
   int num_q = m + 5 * b + 3;
   int offset = 3 * b + 2;
   // int offset = 0;
   int anc = b + offset;
   int x_b = offset;
   int x_1 = anc - 1;
   int x_m = anc - m;
   int w_m_1 = anc + b + m;
   int w_b = anc + 1;
   int w_0 = w_b + b;
   int y_h = x_b;
   int y_0 = y_h - 2; // ym = 2;
   int y_l = y_0 - 2 * b;
   int z_h = y_l;
   int z_0 = z_h; // zm = 0;
   int z_l = z_0 - b;
   Register xr = Register_new(x_b, x_1 + 1, x_1 + 1);
   Register wr = Register_new(w_b, w_0, w_m_1 + 1);
   Register yr = Register_new(y_l, y_0, y_h);
   Register zr = Register_new(z_l, z_0, z_h);
   
   if (VERBOSE) {
      fprintf(stderr, "num_q: %d\n", num_q);
      fprintf(stderr, "anc: %d\n", anc);
      fprintf(stderr, "x_b: %d\n", x_b);
      fprintf(stderr, "x_1: %d\n", x_1);
      fprintf(stderr, "x_m: %d\n", x_m);
      fprintf(stderr, "w_m_1: %d\n", w_m_1);
      fprintf(stderr, "w_b: %d\n", w_b);
      fprintf(stderr, "w_0: %d\n", w_0);
      fprintf(stderr, "int: %d\n", (m - n + 1 + 2 * b));
      // fprintf(stderr, "m: %d\n", m);
   }
   
   State state = State_new_integer(num_q, w << (m - n + 1 + 2 * b + offset));
   printf("\nInput number: %f\t", w / pow(2.0,n-m));
   printf("Actual Inverse: \t%.15f\n-------------------------\n", 
         pow(2.0,n-m) / w);
   State_print(state); fflush(stdout);
   
   state = Inverse_init(state,wr,xr,anc);
   double temp = Register_measure_as_double(xr,state);
   printf("\nInitial approximation: 1/%ld", (long) (1.0 / temp));
   printf("\t\t\t = \t%.15f\n-----------------------------\n", temp);
   State_print(state); fflush(stdout);
   
   
   Operator ms = Operator_multiply_subtract_from(num_q, wr, xr, yr); // y-= wx
   Operator ma = Operator_multiply_add(num_q, xr, yr, zr); // z+= xy
   int s = ceil(log2(b));
   for (int i = 1; i <= s; i++) {
      Register_set(yr, state, 1 << (y_h - y_l - 1)); // y := 2
      state = Operator_apply(ms, state);
      
      // double temp0 = Register_measure_as_double(xr,state);
      // printf("\nx_%d: 1/%ld", i-1, (long) (1.0 / temp0));
      // printf("\t\t\t = \t%.15f\n-----------------------------\n", temp0);
      // State_print(state); fflush(stdout);
      
      if (VERBOSE) {
         double temp1 = Register_measure_as_double(yr,state);
         printf("\n2 - w x_%d: 1 + 1/%.2f", i-1, 
               ((long) (100.0 / (temp1-1)))/100.0);
         printf("\t\t\t = \t%.15f\n-----------------------------\n", temp1);
         State_print(state); fflush(stdout);
      }
      
      Register_set(zr, state, 0L); // z := 0
      state = Operator_apply(ma, state);
      // Register_set(wr, state, 0L); // z := 0
      // Operator ma2 = Operator_multiply_add(num_q, xr, yr, wr); // z+= xy
      // state = Operator_apply(ma2, state);
      
      if (VERBOSE) {
         double temp2 = Register_measure_as_double(zr,state);
         printf("\nx_%d(2 - w x_%d): 1/%.2f", i-1, i-1, 
               ((long) (100.0 / temp2))/100.0);
         printf("\t\t\t = \t%.15f\n-----------------------------\n", temp2);
         State_print(state); fflush(stdout);
      }
      
      Register_set(xr, state, 0L); // x := 0
      for (int j = 0; j < b; j++) { // x := z
         int ctr[2];
         ctr[0] = z_0 - 1 - j;
         ctr[1] = -1;
         int trg[2];
         trg[0] = x_1 - j;
         trg[1] = -1;
         Operator cnot = Operator_cnot(num_q,ctr,trg);
         state = Operator_apply(cnot, state);
         Operator_free(cnot);
      }

      double temp = Register_measure_as_double(xr,state);
      printf("\nx_%d: 1/%.2f", i, ((long) (100.0 / temp))/100.0);
      printf("\t\t\t = \t%.15f\n-----------------------------\n", temp);
      State_print(state); fflush(stdout);
      
   }
   Operator_free(ms);
   Operator_free(ma);
   
   printf("\n");
   return state;
}

void Inverse_run_tests(void) {
   double input;
   printf("Input w.\nw = ");
	scanf("%lf", &input);
   
   int n = 18;
   int m = 8;
   long w = input * pow(2,n-m);
   int b = 10;
   assert(m + 2 * b + 1 < 31); // otherwise, get integer overflow
   State s = Inverse_inverse(w,n,m,b);
}

