/* 
 * vector.h
 * --------
 * A Vector object is a representation of a sparse array/vector of 
 * complex doubles and contains index,value pairs in index order.
 * 
 * Supports many dictionary operations such as inserting elements, deleting 
 * elements, get operation (lookup), mapping across all index,value pairs
 * 
 * Also supports vector operations such as vector sum, 
 * norm (magnitude), scalar product, inner product
 * 
 * Author: Barak Nehoran
 * Advisor: Iasonas Petras
 */


#ifndef VECTOR_H
#define VECTOR_H

#include "pair.h"
#include <stdlib.h>
#include <assert.h>
#include <complex.h>
#include <stdbool.h>

/* A Vector object is a sparse array/vector of complex floats */

typedef struct VectorStruct* Vector;

/*--------------------------------------------------------------------*/

/* Return an empty Vector object. */

Vector Vector_new_empty(void);

/* Return a new Vector object, or NULL if insufficient memory is
   available. */

Vector Vector_new(long index, double complex val);

/* Variant of Vector_new() with a bignum index */
Vector Vector_new_bn(const bignum index, double complex val);

/*--------------------------------------------------------------------*/

/* Free vec. */

void Vector_free(Vector vec);

/*--------------------------------------------------------------------*/

/* Insert val at position index. Return the new vector if successful, 
   or NULL if insufficient memory is available. */

Vector Vector_insert(Vector vec, long index, double complex val);

/* Variant of Vector_insert() with a bignum index */
Vector Vector_insert_bn(Vector vec, const bignum index, double complex val);

/*--------------------------------------------------------------------*/

/* Adds val at position index. Return the new vector if successful, 
   or NULL if insufficient memory is available. */

Vector Vector_addto(Vector vec, long index, double complex val);

/* Variant of Vector_addto() with a bignum index */
Vector Vector_addto_bn(Vector vec, const bignum index, double complex val);

/*--------------------------------------------------------------------*/

/* Returns the value of vec at the specified index */

double complex Vector_get(const Vector vec, long index);

/* Variant of Vector_get() with a bignum index */
double complex Vector_get_bn(const Vector vec, const bignum index);

/*--------------------------------------------------------------------*/

/* Deletes the value at the given index. Returns the new vector. */

Vector Vector_delete(Vector vec, long index);

/* Variant of Vector_delete() with a bignum index */
Vector Vector_delete_bn(Vector vec, const bignum index);

/*--------------------------------------------------------------------*/

/* Return true if vec is empty, or false otherwise. */

bool Vector_isEmpty(const Vector vec);

/*--------------------------------------------------------------------*/

/* Return the number of explicitly represented (mostly nonzero) elements 
   in the vector. Includes elements that are explicitly represented 
   but may have value zero. 
   WARNING: the non-bignum version can have overflow if the result is
   larger than can be stored in a long. */

long Vector_nnz(const Vector vec);

/* A safer variant of Vector_nnz() that stores the result as a bignum 
   in the argument named result. */
void Vector_nnz_bn(bignum result, const Vector vec);

/*--------------------------------------------------------------------*/

/* Return an array containing the index,value pairs in the vector
   in ascending order. The array is terminated with a -1 index.
   Note: The array returned placed in the heap and should be freed 
   by the client.
   Note: May return NULL if there are too many elements to store in 
   an array.
   
   WARNING: Vulnerable. Should only be used for testing! 
   Pairs in the returned array should not be modified by the client
   as they are they could invalidate the vector and cause undefined
   behavior. */

// struct Pair* Vector_contents(Vector vec);

/*--------------------------------------------------------------------*/

/* Apply f to every component of the vector and return the result as a
   new vector. Return NULL if insufficient memory available. */

Vector Vector_map(Vector vec, struct Pair (*f)(struct Pair, void*),
        void* extra);

/* Same as map, but where the function may map one component to multiple.
   The Pair array returned by f must be terminated by a -1 index. */

Vector Vector_map_mult(Vector vec, struct Pair* (*f)(struct Pair, void*), 
        void* extra);

/*--------------------------------------------------------------------*/

// /* Deprecated: Return a reduced version of the vector with all 
//    zero and almost-zero values removed. */

// Vector Vector_reduce(Vector vec);

/*--------------------------------------------------------------------*/

/* Return the vector with all values multiplied by the scalar c. */

Vector Vector_scalar_prod(Vector vec, double complex c);

/*--------------------------------------------------------------------*/

/* Return a vector with values that are the complex conjugates of 
   those of vec. */

Vector Vector_conj(Vector vec);

/*--------------------------------------------------------------------*/

/* Return a vector that is the sum of the two vectors. 
   Does not modify either argument vector. */

Vector Vector_sum(const Vector v1, const Vector v2);

/*--------------------------------------------------------------------*/

/* Return a vector that is the element-wise product of the two vectors. 
   Does not modify either argument vector. */

Vector Vector_element_prod(const Vector v1, const Vector v2);

/*--------------------------------------------------------------------*/

/* Return the inner product of the two vectors. 
   Does not modify either argument vector. */

double complex Vector_inner_prod(const Vector v1, const Vector v2);

/*--------------------------------------------------------------------*/

/* Return the L2 norm or magnitude of a vector.
   Does not modify the argument vector. */

double complex Vector_norm(const Vector vec);

/*--------------------------------------------------------------------*/

/* Prints out the Pairs contained in the vector, one per line. */

void Vector_print(const Vector vec);

/*--------------------------------------------------------------------*/

/* Runs the unit testing code. Except when preblems arise, running 
   this function should produce no visible effect. */
   
void Vector_run_tests (int n);

/*--------------------------------------------------------------------*/


#endif