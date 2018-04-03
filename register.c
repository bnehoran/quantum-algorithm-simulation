/* 
 * register.c
 * ----------
 * A Register is a definition of a string of qubits as designating 
 * a fixed point binary number.
 * 
 * Registers are represented by a trio of indices low,zero,high 
 * representing the indices of the qubits representing the register.
 * low is the lowest index of the register.
 * zero is the index of the qubit representing the 2^0 bit.
 * high is the lowest index just above the register. That is, 
 * the lowest index of the next register.
 * 
 * This is defined such that high - low is the total number of 
 * qubits in the register, high - zero is the number of 
 * qubits in the integer part of the representation, and 
 * zero - low is the number of qubits in the fractional part.
 * 
 * Author: Barak Nehoran
 * Advisor: Iasonas Petras
 */

#include "register.h"
#include "state.h"
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

// struct RegisterStruct {
//     int low;
//     int zero;
//     int high;
// };

/* Creates and returns a new register with the given low, zero, and high
   indices. \n
   Requirement: low <= zero <= high */
Register Register_new(int low, int zero, int high) {
   Register reg = (Register)malloc(sizeof(struct RegisterStruct));
   memset(reg, 0x0, sizeof(struct RegisterStruct));
   reg->low = low;
   reg->zero = zero;
   reg->high = high;
   return reg;
}

/* Frees the memory held by reg. */
void Register_free(Register reg) {
   free(reg);
}

/* Returns whether or not the two registers r1 and r2 intersect.
   Intersecting registers can be an indication of potential problems. */
bool Register_have_no_overlap(Register r1, Register r2) {
   assert(Register_isValid(r1) && Register_isValid(r2));
   return (r1->low >= r2->high || r2->low >= r1->high);
}

/* Initializes the register within the state s to the binary 
   representation of the integer x. */
State Register_set(Register reg, State s, long x) {
   return State_set_register (s, x, reg->low, reg->high);
}

/* Measures the value of the register reg in the state s and 
   returns it either as an integer (including both the integer
   and fractional parts of the register) or as a floating point
   number that represents the actual value (but which may or 
   may not have roundoff error in the fractional part). */
long Register_measure(Register reg, State s) {
   long x = State_measure_all(s);
   return (x % (1L << reg->high)) >> reg->low;
}

double Register_measure_as_double(Register reg, State s) {
   long x = Register_measure(reg,s);
   return x / pow(2.0, reg->zero - reg->low);
}

/* Returns whether the register reg is valid. */
bool Register_isValid(Register reg) {
   return (0 <= reg->low && reg->low <= reg->zero && reg->zero <= reg->high);
}

/*--------------------------------------------------------------------*/
/* Unit testing for register.c                                        */
/*--------------------------------------------------------------------*/

void test_constructor (int n) {}
void test_isValid (int n) {}
void test_have_no_overlap (int n) {}
void test_set (int n) {}
void test_measure (int n) {}
void test_measure_as_double (int n) {}

/* Runs the unit testing code. Except when preblems arise, running 
   this function should produce no visible effect. */
void Register_run_tests(int n) {
   test_constructor (n);
   test_isValid (n);
   test_have_no_overlap (n);
   test_set (n);
   test_measure (n);
   test_measure_as_double (n);
}
