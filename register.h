/* 
 * register.h
 * ----------
 * A Register is a definition of a string of qubits as designating 
 * a fixed precision binary number.
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

#ifndef REGISTER_H
#define REGISTER_H

#include "state.h"
#include "pair.h"
#include <stdbool.h>

/* Registers are represented by a trio of indices low,zero,high representing
   the indices of the qubits representing the register.
   low is the lowest index of the register.
   zero is the index of the qubit representing the 2^0 bit.
   high is the lowest index just above the register. That is, 
   the lowest index of the next register. */
typedef struct RegisterStruct* Register;

struct RegisterStruct {
    int low;
    int zero;
    int high;
};

/* Creates and returns a new register with the given low, zero, and high
   indices. \n
   Requirement: low <= zero <= high */
Register Register_new(int low, int zero, int high);

/* Frees the memory held by reg. */
void Register_free(Register reg);

/* Returns whether or not the two registers r1 and r2 intersect.
   Intersecting registers can be an indication of potential problems. */
bool Register_have_no_overlap(Register r1, Register r2);

/* Initializes the register within the state s to the binary 
   representation of the integer x. */
State Register_set(Register reg, State s, long x);

/* Measures the value of the register reg in the state s and 
   returns it either as an integer (including both the integer
   and fractional parts of the register) or as a floating point
   number that represents the actual value (but which may or 
   may not have roundoff error in the fractional part). */
long Register_measure(Register reg, State s);
double Register_measure_as_double(Register reg, State s);

/* Returns whether the register reg is valid. */
bool Register_isValid(Register reg);

/* Runs the unit testing code. Except when preblems arise, running 
   this function should produce no visible effect. */
void Register_run_tests(int n);

#endif