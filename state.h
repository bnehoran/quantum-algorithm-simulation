/* 
 * state.h
 * --------
 * A State object is a vector reprsentation of the quantum mechanical state
 * of an n-qubit quantum computer.
 * 
 * Supports adding and removing qubits, as well as several measurement 
 * operations, and mapping operations to convert from one State to another.
 * 
 * A state should always be normalized.
 * 
 * Author: Barak Nehoran
 * Advisor: Iasonas Petras
 */

#ifndef STATE_H
#define STATE_H

#include "pair.h"
#include <stdlib.h>
#include <assert.h>
#include <complex.h>
#include <stdbool.h>

typedef struct StateStruct* State;

enum Mode {QUANTUM, COMPUT_BASIS};

/*
 * Standard Useful Values / Constructors
 */


/* Returns a basis State that is a n - qubit
   representation of the integer.
   integer must be at least 0 and at most 2^n - 1 */
State State_new_integer (int n, long integer);

/* Gives a State of n qubits all initialized to |0>. */
State State_new_zero (int n);

/* Returns the single-qubit State representing a|0> + b|1> */
State State_new_single_qubit (double complex a, double complex b);

/* Returns the a positive equal superposition of all possible n-qubit states. */
State State_new_hadamard (int n);

/* Frees the state s. */
void State_free (State s);

/* Return the number of cubits represented by the state. */
int State_num_qubits (State s);

/*
 * Modifying a state
 */

/* Takes a State |s> of n qubits and adds a qubit to the State
   at position pos <= n to give a State of n+1 qubits
   The new qubit is initialized in the State
   a|0> + b|1>, or |0> by default.
   The default position is at the end of the list of qubits. */
State State_add_qubit (State s);
State State_add_qubit_pos (State s, int pos);
State State_add_qubit_as (State s, double complex a, double complex b, int pos);
// State State_add_qubit (State s; [double complex a; double complex b;]);

/* Takes a State |s> of n qubits and removes the qubit
   at position pos < n to give a State of n-1 qubits.
   This qubit must be completely factorable from the rest of the State
   to prevent the situation in which classical probabilities arise. */
State State_forget_qubit (State s, int pos);

/* Returns whether the qubit in position pos can be
   factored out as a tensor product.
   Additionally, as_a_zero specifies that this qubit should
   factor out as the State |0> */
bool State_is_factorable (State s, int pos);
bool State_is_factorable_as_zero (State s, int pos);

/* returns the inner product <phi|psi> of two States |phi> and |psi> */
double State_inner_product (State phi, State psi);

/* Subroutine to give the L2 norm of the State |s>
   Should always return 1.0 for a normalized State */
double State_norm (State s);

/* Returns a normalized (norm(s) == 1.0) version of s */
State State_normalize (State s);

/* Takes a State |phi> of n qubits and a State |psi> of m qubits
   and returns the State |phi>|psi> of n+m qubits */
State State_tensor_product (State phi, State psi);

/* Apply f to every component of the state and return the result as a
   new state. Return NULL if insufficient memory available. */
State State_map (State s, struct Pair (*f)(struct Pair, void*), void* extra);

/* Same as map, but where the function may map one component to multiple.
   The Pair array returned by f must be terminated by a -1 index. */
State State_map_mult (State s, struct Pair* (*f)(struct Pair, void*), 
      void* extra);

//
// Measurement
//

/* Takes an n-qubit State |s> and collapses it to the State
   produced after measuring mes on the k-th qubit.
   Returns the resulting n - qubit State.
   The parameter mes must be 0 or 1 */
State State_collapse (State s, int k, int mes);


/* Takes an n-qubit State |s> and measures all n qubits
   Returns the result as an n-bit integer */
long State_measure_all (State s);


/* Takes an n-qubit State |s> and "measures" the k-th qubit
   Returns an 0|1 integer representing the measurement result */
long State_measure (State s, int k);

/* Initializes a k-qubit register within state s to the binary 
   representation of the integer x. The register begins at qubit 
   number x_0 and ends at qubit number x_k-1. (x_k is the first
   qubit of the next register, or the number of qubits if the 
   register extends to the end of the state. */
State State_set_register (State s, long x, int x_0, int x_k);

/* Prints out the state in the computational basis with 
   basis state, amplitude Pairs printed one per line.
   Mainly to be used for testing purposes. */
void State_print(State s);

/* Runs the unit testing code. Except when preblems arise, running 
   this function should produce no visible effect. */
void State_run_tests(int n);



#endif