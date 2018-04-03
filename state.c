/* 
 * state.c
 * --------
 * A State object is a vector reprsentation of the quantum mechanical state
 * of an n-qubit quantum computer.
 * It uses a Vector (vector.h) as a base for representing a state
 * of n qubits.
 * 
 * Supports adding and removing qubits, as well as several measurement 
 * operations, and mapping operations to convert from one State to another.
 * 
 * A state should always be normalized.
 * 
 * Author: Barak Nehoran
 * Advisor: Iasonas Petras
 */

#include "state.h"
#include "vector.h"
#include <stdlib.h>
#include <assert.h>
#include <complex.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

struct StateStruct {
   // the vector representing the quantum state
   Vector vector;
   // then number of quibits represented by this state
   int n;
   // can be either QUANTUM or COMPUT_BASIS representing whether
   // superpositions of basis states are allowed or not
   enum Mode mode;
};


/*
 * Standard Useful Values / Constructors
 */

/* Returns a basis State that is a n - qubit
   representation of the integer.
   integer must be at least 0 and at most 2^n - 1 */
State State_new_integer (int n, long integer) {
   State s = (State)malloc(sizeof(struct StateStruct));
   if (s == NULL) // insufficient memory
      return NULL;
   s->n = n;
   s->vector = Vector_new(integer,1.0);
   s->mode = COMPUT_BASIS;
   return s;
}

State State_new_integer_bn (int n, bignum integer) {
   State s = (State)malloc(sizeof(struct StateStruct));
   if (s == NULL) // insufficient memory
      return NULL;
   s->n = n;
   s->vector = Vector_new_bn(integer,1.0);
   s->mode = COMPUT_BASIS;
   return s;
}

/* Gives a State of n qubits all initialized to |0>. */
State State_new_zero (int n) {
    return State_new_integer(n,0);
}

/* Returns the single-qubit State representing a|0> + b|1> */
State State_new_single_qubit (double complex a, double complex b) {
   State s = (State)malloc(sizeof(struct StateStruct));
   if (s == NULL) // insufficient memory
      return NULL;
   double complex mag = csqrt(conj(a)*a + conj(b)*b);
   a /= mag; // ensure that the new state is normalized.
   b /= mag;
   s->n = 1;
   s->vector = Vector_new_empty();
   s->vector = Vector_insert(s->vector, 0, a);
   s->vector = Vector_insert(s->vector, 1, b);
   s->mode = QUANTUM;
   return s;
}

/* Returns the a positive equal superposition of all possible n-qubit states.
   Uses exponential time and space. Not recommended for large n. */

State State_new_hadamard (int n) {
   State s = (State)malloc(sizeof(struct StateStruct));
   if (s == NULL) // insufficient memory
      return NULL;
   s->n = n;
   s->vector = Vector_new_empty();
   long N = 1 << n;
   assert(N > 1);
   double complex amplitude = 1.0 / sqrt(N);
   for (long i = 0; i < N; i++) {
      s->vector = Vector_insert(s->vector, i, amplitude);
   }
   s->mode = QUANTUM;
   return s; // state is already normalized, so judt return it.
}

/* Frees the state s. */
void State_free (State s) {
   Vector_free(s->vector);
   free(s);
}

/* Return the number of cubits represented by the state. */
int State_num_qubits (State s) {
   return s->n;
}

/*
 * Modifying a state
 */

/* Takes a State |s> of n qubits and adds a qubit to the State
   at position pos <= n to give a State of n+1 qubits
   The new qubit is initialized in the State
   a|0> + b|1>, or |0> by default.
   The default position is at the end of the list of qubits. */
State State_add_qubit (State s) {
   long kth_bit = 1;
   struct Pair inesrt_bit(struct Pair p, void* extra) {
      p.index = p.index << 1;
      return p;
   }
   Vector temp = s->vector;
   s->vector = Vector_map(s->vector, &inesrt_bit, NULL);
   free(temp);
   s->n = s->n + 1;
   return s;
}

// TODO: Expand so that able to add in pos >= 63
State State_add_qubit_pos (State s, int pos) {
   assert(pos < 63);
   long kth_bit  = 1 << pos;
   struct Pair inesrt_bit(struct Pair p, void* extra) {
      long i = p.index;
      long j_bot = i % kth_bit;
      long j_top = i - j_bot;
      long j = (j_top << 1) + j_bot;
      p.index = j;
      return p;
   }
   Vector temp = s->vector;
   s->vector = Vector_map(s->vector, &inesrt_bit, NULL);
   free(temp);
   s->n = s->n + 1;
   return s;
}

// TODO: Expand so that able to add in pos >= 63
State State_add_qubit_as (State s, double complex a, double complex b, int pos){
   assert(pos < 63);
   double complex mag = csqrt(conj(a)*a + conj(b)*b);
   a /= mag; // ensure that the new state is normalized.
   b /= mag;
   long kth_bit  = 1 << pos;
   struct Pair arr[3];
   struct Pair* inesrt_bit(struct Pair p, void* extra) {
      long i = p.index;
      long j_bot = i % kth_bit;
      long j_top = i - j_bot;
      long j = (j_top << 1) + j_bot;
      arr[0].index = j;  // with a 0 in the kth bith
      arr[0].val = a * p.val;
      arr[1].index = j + kth_bit; // with a 1 in the kth bit
      arr[1].val = b * p.val;
      arr[2].index = -1;
      arr[2].val = 0.0;
      return arr;
   }
   Vector temp = s->vector;
   s->vector = Vector_map_mult(s->vector, &inesrt_bit, NULL);
   free(temp);
   s->n = s->n + 1;
   s->mode = QUANTUM;
   return s;
}

// State State_add_qubit (State s; [double complex a; double complex b;]);

/* Takes a State |s> of n qubits and removes the qubit
   at position pos < n to give a State of n-1 qubits.
   This qubit must be completely factorable from the rest of the State
   to prevent the situation in which classical probabilities arise. */
State State_forget_qubit (State s, int pos) {
   // TODO
   return NULL;
}

/* Returns whether the qubit in position pos can be
   factored out as a tensor product.
   Additionally, as_a_zero specifies that this qubit should
   factor out as the State |0> */
bool State_is_factorable (State s, int pos) {
   //TODO
   return false;
}

bool State_is_factorable_as_zero (State s, int pos) {
   //TODO
   return false;
}

/* returns the inner product <phi|psi> of two States |phi> and |psi> */
double State_inner_product (State phi, State psi) {
   return Vector_inner_prod(phi->vector, psi->vector);
}

/* Subroutine to give the L2 norm of the State |s>
   Should always return 1.0 for a normalized State */
double State_norm (State s) {
   return Vector_norm(s->vector);
}

/* Returns a normalized (norm(s) == 1.0) version of s */
State State_normalize (State s) {
   double factor = Vector_norm(s->vector);
   s->vector =  Vector_scalar_prod(s->vector, 1.0 / factor);
   return s;
}

/* Takes a State |phi> of n qubits and a State |psi> of m qubits
   and returns the State |phi>|psi> of n+m qubits */
State State_tensor_product (State phi, State psi) {
   //TODO
   return NULL;
}

/* Apply f to every component of the state and return the result as a
   new state. Return NULL if insufficient memory available. */
State State_map (State s, struct Pair (*f)(struct Pair, void*), void* extra) {
   Vector v_new = Vector_map(s->vector,f,extra);
   if (v_new == NULL) // insufficient memory
      return NULL;
   State s_new = (State)malloc(sizeof(struct StateStruct));
   if (s_new == NULL) // insufficient memory
      return NULL;
   s_new->n = s->n;
   s_new->vector = v_new;
   s_new->mode = s->mode;
   return s_new;
}

/* Same as map, but where the function may map one component to multiple.
   The Pair array returned by f must be terminated by a -1 index. */
State State_map_mult (State s, struct Pair* (*f)(struct Pair, void*), 
      void* extra) {
   assert (s->vector != NULL); // invalid state
   Vector v_new = Vector_map_mult(s->vector,f,extra);
   if (v_new == NULL) // insufficient memory OR empty vector
      return NULL;
   State s_new = (State)malloc(sizeof(struct StateStruct));
   if (s_new == NULL) // insufficient memory
      return NULL;
   s_new->n = s->n;
   s_new->vector = v_new;
   s_new->mode = QUANTUM;
   return s_new;
}

//
// Measurement
//

/* Takes an n-qubit State |s> and collapses it to the State
   produced after measuring mes on the k-th qubit.
   Returns the resulting n - qubit State.
   The parameter mes must be 0 or 1 */
State State_collapse (State s, int k, int mes) {
   assert(mes == 0 || mes == 1);
   // Vector result = Vector_new_empty();
   long kth_bit = 1 << k;
   struct Pair coll(struct Pair p, void* extra) {
      if ((mes && (p.index & kth_bit)) || (!mes && !(p.index & kth_bit))) {
         return p;
      }
      p.index = -1;
      p.val = 0.0;
      return p;
   }
   Vector temp = s->vector;
   s->vector = Vector_map(temp, &coll, NULL);
   Vector_free(temp);
   return State_normalize(s);
}


/* Takes an n-qubit State |s> and measures all n qubits
   Returns the result as an n-bit integer */
long State_measure_all (State s) {
   double r = rand() / (RAND_MAX + 0.0);
   double* r_pointer = &r;
   long result = -1;
   long* result_pointer = &result; 
   struct Pair measure_helper(struct Pair p, void* extra) {
      double prob = creal(conj(p.val)*p.val);
      *(r_pointer) -= prob;
      if (*(r_pointer) <= 0 && *(result_pointer) == -1) {
         *(result_pointer) = p.index;
      }
      p.index = -1;
      return p;
   }
   Vector_free(Vector_map(s->vector, &measure_helper, NULL));
   return result;
}


/* Takes an n-qubit State |s> and "measures" the k-th qubit
   Returns an 0|1 integer representing the measurement result */
long State_measure (State s, int k) {
   long kth_bit = 1 << k;
   return (State_measure_all(s) & kth_bit);
}

/* Initializes a k-qubit register within state s to the binary 
   representation of the integer x. The register begins at qubit 
   number x_0 and ends at qubit number x_k-1. (x_k is the first
   qubit of the next register, or the number of qubits if the 
   register extends to the end of the state.) */
State State_set_register (State s, long x, int x_0, int x_k) {
   assert(x >= 0);
   assert(x_k >= x_0);
   assert(x >> (x_k - x_0) == 0L); // x should fit in the register
   if (x_k == x_0) return s;
   // Vector result = Vector_new_empty();
   // long kth_bit = 1 << k;
   struct Pair set_register_helper(struct Pair p, void* extra) {
      long prev_contents = (p.index % (1L << x_k)) >> x_0;
      p.index += (x - prev_contents) << x_0;
      return p;
   }
   Vector temp = s->vector;
   s->vector = Vector_map(temp, &set_register_helper, NULL);
   Vector_free(temp);
   return s;
}

State State_set_register_bn (State s, long x, int x_0, int x_k) {
   assert(x >= 0);
   assert(x_k >= x_0);
   assert(x >> (x_k - x_0) == 0L); // x should fit in the register
   if (x_k == x_0) return s;
   // Vector result = Vector_new_empty();
   // long kth_bit = 1 << k;
   struct Pair set_register_helper(struct Pair p, void* extra) {
      long prev_contents = (p.index % (1L << x_k)) >> x_0;
      p.index += (x - prev_contents) << x_0;
      return p;
   }
   Vector temp = s->vector;
   s->vector = Vector_map(temp, &set_register_helper, NULL);
   Vector_free(temp);
   return s;
}

/* Prints out the state in the computational basis with 
   basis state, amplitude Pairs printed one per line.
   Mainly to be used for testing purposes. */
void State_print(State s) {
   printf("State of %d qubits ",s->n);
   if (s->mode == COMPUT_BASIS) {
      printf("in the basis state:\n");
   }
   else {
      printf("in a superposition of %ld basis states:\n",Vector_nnz(s->vector));
   }
   // Vector_print(s->vector); // alternative to just print as a vector
   struct Pair print_pair (struct Pair p, void* extra) {
      char bits[s->n + 1];
      bits[s->n] = 0x0;
      for (int i = 0; i < s->n; i++) {
         if (p.index & (1L << (s->n - 1 - i)))
            bits[i] = '1';
         else
            bits[i] = '0';
      }
      printf("(%.1f+%.1fi) |%s (%lx)>\n", 
            creal(p.val), cimag(p.val), bits, p.index);
      p.index = -1;
      return p;
   }
   
   Vector_free(Vector_map(s->vector,&print_pair,NULL)); 
   printf("\n");
   fflush(stdout);
}


/*--------------------------------------------------------------------*/
/* Unit testing for state.c                                          */
/*--------------------------------------------------------------------*/

//TODO: write tests (see tests in vector.c)
bool State_isValid(State s){
   //TODO
   return false;
}

void test_State_new_integer (int n) {
   
   // State_new_integer (int n, long integer);
}
void test_State_new_zero (int n) {
   
   // State_new_zero (int n);
}
void test_State_new_single_qubit (int n) {
   
   // State_new_single_qubit (double complex a, double complex b);
}
void test_State_add_qubit (int n) {
   
   // State_add_qubit (State s);
}
void test_State_forget_qubit (int n) {
   
   // State_forget_qubit (State s, int pos);
}
void test_State_is_factorable (int n) {
   
   // State_is_factorable (State s, int pos);
}
void test_State_inner_product (int n) {
   
   // State_inner_product (State phi, State psi);
}
void test_State_norm (int n) {
   
   // State_norm (State s);
}
void test_State_normalize (int n) {
   
   // State_normalize (State s);
}
void test_State_tensor_product (int n) {
   
   // State_tensor_product (State phi, State psi);
}
void test_State_collapse (int n) {
   
   // State_collapse (State s, int k, int mes);
}
void test_State_measure_all (int n) {
   
   // State_measure_all (NULL);
}
void test_State_measure (int n) {
   
   // State_measure (State s, int k);
}

/* Runs the unit testing code. Except when preblems arise, running 
   this function should produce no visible effect. */
void State_run_tests(int n) {
   test_State_new_integer (n);
   test_State_new_zero (n);
   test_State_new_single_qubit (n);
   test_State_add_qubit (n);
   test_State_forget_qubit (n);
   test_State_is_factorable (n);
   test_State_inner_product (n);
   test_State_norm (n);
   test_State_normalize (n);
   test_State_tensor_product (n);
   test_State_collapse (n);
   test_State_measure_all (n);
   test_State_measure (n);
}

int State_main (void) {
   srand(time(NULL));

   int n = 10; // Number of qubits to use for tests
   State_run_tests(n);
   
   // TODO finish transition to bignum
   
   // int n = 8;
   // State s = State_new_hadamard(n);
   // struct Pair print_pair (struct Pair p, void* extra) {
   //    printf("%X: %.1f+%.1fi\n", 
   //          p.index, creal(p.val), cimag(p.val));
   //    p.index = -1;
   //    return p;
   // }
   // // Vector_map(s->vector, &print_pair, NULL);
   // int res = State_measure_all(s);
   // printf("measured: %d out of 2^%d = %d\n", res, n, 1<<n);
   // for (int i = 0; i < n; i++) {
   //    int mes = (State_measure(s,i) && true);
   //    printf("measured: %dth qubit and got a %d\n", i, mes);
   //    s = State_collapse(s,i,mes);
   //    assert(isValid(s->vector));
   //    Vector_map(s->vector, &print_pair, NULL);
   // }
   // printf("collapse!!!!\n");
   // s = State_collapse(s,1,0);
   // Vector_map(s->vector, &print_pair, NULL);
   
   // int i = p.index;
   // int j_bot = i % kth_bit;
   // int j_top = i - mes * kth_bit - j_bot;
   // int j = j_top >> 1 + j_bot;
   // p.index = j
   
   return 0;
}