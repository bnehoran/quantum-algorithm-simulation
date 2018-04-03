/* 
 * vector.c
 * --------
 * A Vector object is a sparse array/vector of complex doubles
 * represented as a linked list containing index,value pairs in index order.
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

#include "vector.h"
#include <stdlib.h>
#include <assert.h>
#include <complex.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

/* VectorStruct is the concrete representation of the Vector type.
   It also doubles as a node in the linked list. 
   An empty vector is represented by a NULL pointer. (This is, of 
   course, abstracted from the client, however, so the client should 
   not assume this behavior.) */

struct VectorStruct {
   struct Pair p;
   struct VectorStruct* next;
};

/*--------------------------------------------------------------------*/

/* Return an empty Vector object.
   Note: an empty vector is represented by a NULL pointer. */

Vector Vector_new_empty(void) {
   return NULL;
}

/* Return a new Vector object, or NULL if insufficient memory is
   available. */

Vector Vector_new(long index, double complex val) {
   bn_new_with(_n, index);
   Vector result = Vector_new_bn(_n, val);
   bn_clear(_n);
   return result;
}

/* Variant of Vector_new() with a bignum index */
Vector Vector_new_bn(const bignum index, double complex val) {
   assert(bn_sgn(index) >= 0); // only positive indices allowed
   Vector vec = (Vector)malloc(sizeof(struct VectorStruct));
   if (vec == NULL) return NULL;
   bn_initialize(vec->p.index_bn);
   bn_set(vec->p.index_bn, index); // vec->p.index = index
   vec->p.index = bn_lval(index);
   vec->p.val = val;
   vec->next = NULL;
   return vec;
}

/*--------------------------------------------------------------------*/

/* Free vec. */

void Vector_free(Vector vec) {
   // assert(vec); // remove: vec == NULL for empty vector
   Vector cur = vec;
   Vector nex = NULL;
   while (cur != NULL) {
      nex = cur->next;
      bn_clear(cur->p.index_bn);
      free(cur);
      cur = nex;
   }
}

/*--------------------------------------------------------------------*/

/* Insert val at position index. Return the new vector if successful, 
   or NULL if insufficient memory is available. */

Vector Vector_insert(Vector vec, long index, double complex val) {
   bn_new_with(_n, index);
   Vector result = Vector_insert_bn(vec, _n, val);
   bn_clear(_n);
   return result;
}

/* Variant of Vector_insert() with a bignum index */
Vector Vector_insert_bn(Vector vec, const bignum index, double complex val) {
   if (Vector_isEmpty(vec)) return Vector_new_bn(index,val); // empty case
   Vector cur = vec;
   Vector prev = NULL;
   // for (cur = vec; cur != NULL && cur->p.index < index; 
   for (cur = vec; cur != NULL && bn_cmp(cur->p.index_bn, index) < 0; 
         prev = cur, cur = cur->next);
   // if (cur != NULL && cur->p.index == index) { // found index
   if (cur != NULL && bn_cmp(cur->p.index_bn, index) == 0) { // found index
      // replace the existing value with val
      cur->p.val = val;
      return vec;
   }
   // not found, so insert a new node
   Vector new_node = Vector_new_bn(index,val);
   if (new_node == NULL) return NULL; // insufficient memory
   if (cur == NULL) { // insert at the end
      assert(prev != NULL); // cannot both be NULL if vec is nonempty
      prev->next = new_node;
   } 
   else if (prev == NULL) { // insert as first node
      new_node->next = cur;
      vec = new_node;
   }
   else { // insert as interior node
      new_node->next = cur;
      prev->next = new_node;
   }
   return vec;
}

/*--------------------------------------------------------------------*/

/* Adds val at position index. Return the new vector if successful, 
   or NULL if insufficient memory is available. */

Vector Vector_addto(Vector vec, long index, double complex val) {
   bn_new_with(_n, index);
   Vector result = Vector_addto_bn(vec, _n, val);
   bn_clear(_n);
   return result;
}

/* Variant of Vector_addto() with a bignum index */
Vector Vector_addto_bn(Vector vec, const bignum index, double complex val) {
   if (Vector_isEmpty(vec)) return Vector_new_bn(index,val); // empty case
   Vector cur = vec;
   Vector prev = NULL;
   // for (cur = vec; cur != NULL && cur->p.index < index; 
   for (cur = vec; cur != NULL && bn_cmp(cur->p.index_bn, index) < 0; 
         prev = cur, cur = cur->next);
   // if (cur != NULL && cur->p.index == index) { // found index
   if (cur != NULL && bn_cmp(cur->p.index_bn, index) == 0) { // found index
      // add val to the existing value
      cur->p.val = cur->p.val + val; //only line different than insert
      return vec;
   }
   // not found, so insert a new node
   Vector new_node = Vector_new_bn(index,val);
   if (new_node == NULL) return NULL; // insufficient memory
   if (cur == NULL) { // insert at the end
      assert(prev != NULL); // cannot both be NULL if vec is nonempty
      prev->next = new_node;
   } 
   else if (prev == NULL) { // insert as first node
      new_node->next = cur;
      vec = new_node;
   }
   else { // insert as interior node
      new_node->next = cur;
      prev->next = new_node;
   }
   return vec;
}

// TODO: consider factoring out common code from insert() and addto().

/*--------------------------------------------------------------------*/

/* Returns the value of vec at the specified index. Indices that are not 
   explicitly specified are assumed to be zero. */

double complex Vector_get(const Vector vec, long index) {
   bn_new_with(_n, index);
   double complex result = Vector_get_bn(vec, _n);
   bn_clear(_n);
   return result;
}

/* Variant of Vector_get() with a bignum index */
double complex Vector_get_bn(const Vector vec, const bignum index) {
   Vector cur;
   // for (cur = vec; cur != NULL && cur->p.index < index; cur = cur->next);
   for (cur = vec; cur != NULL && bn_cmp(cur->p.index_bn, index) < 0; 
         cur = cur->next);
   // if (cur != NULL && cur->p.index == index) { // found index
   if (cur != NULL && bn_cmp(cur->p.index_bn, index) == 0) { // found index
      return cur->p.val;
   }
   return 0.0; // not found. Default value is 0.0
}

/*--------------------------------------------------------------------*/

/* Deletes the value at the given index. Returns the new vector. */

Vector Vector_delete(Vector vec, long index) {
   bn_new_with(_n, index);
   Vector result = Vector_delete_bn(vec, _n);
   bn_clear(_n);
   return result;
}

/* Variant of Vector_delete() with a bignum index */
Vector Vector_delete_bn(Vector vec, const bignum index) {
   if (Vector_isEmpty(vec)) return vec; // empty case
   Vector cur = vec;
   Vector prev = NULL;
   // for (cur = vec; cur != NULL && cur->p.index < index; 
   for (cur = vec; cur != NULL && bn_cmp(cur->p.index_bn, index) < 0; 
         prev = cur, cur = cur->next);
   // if (cur != NULL && cur->p.index == index) { // found index
   if (cur != NULL && bn_cmp(cur->p.index_bn, index) == 0) { // found index
      if (prev == NULL) { // first node
         vec = cur->next;
         cur->next = NULL;
         Vector_free(cur);
      }
      else { // interior node
         prev->next = cur->next;
         cur->next = NULL;
         Vector_free(cur);
      }
   }
   return vec;
}

/*--------------------------------------------------------------------*/

/* Return true if vec is empty, or false otherwise. */

bool Vector_isEmpty(const Vector vec) {
   return (vec == NULL);
}

/*--------------------------------------------------------------------*/

// // For Reference. Included already in vector.h
// struct Pair {
//    long index;
//    bignum index_bn;
//    double complex val;
// };

/*--------------------------------------------------------------------*/

/* Return the number of explicitly represented (mostly nonzero) elements 
   in the vector. Includes elements that are explicitly represented 
   but may have value zero. 
   WARNING: the version returning a long can have overflow if the result is
   larger than can be stored in a long. If this is a possibility, consider
   using the _bn or _d versions to have the number returned as a bignum
   or as a double, respectively. */

long Vector_nnz(Vector vec) {
   bn_new(_n);
   Vector_nnz_bn(_n, vec);
   long result = bn_lval(_n);
   bn_clear(_n);
   return result;
}

/* A safer variant of Vector_nnz() that stores the result as a bignum 
   in the argument named result. */
void Vector_nnz_bn(bignum result, const Vector vec) {
   bn_setl(result, 0); // result = 0;
   for (Vector cur = vec; cur != NULL; cur = cur->next) {
      bn_incr(result); // result++;
   }
}

/* Safe variant of Vector_nnz_bn() returning a double.
   Some precision may be lost when using this version as it drops any less
   significant bits that will not fit on a double. */
double Vector_nnz_d(Vector vec) {
   bn_new(_n);
   Vector_nnz_bn(_n, vec);
   double result = bn_dval(_n);
   bn_clear(_n);
   return result;
}

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
struct Pair* Vector_contents(Vector vec) {
   bn_new(n); // bignum n;
   Vector_nnz_bn(n, vec); // n = Vector_nnz(vec);
   if (!bn_fits_in_long(n)) return NULL; // too many elements
   long nn = bn_lval(n);
   struct Pair* arr = (struct Pair*)malloc((nn + 1) * sizeof(struct Pair));
   long i = 0;
   for (Vector cur = vec; cur != NULL; cur = cur->next, i++) {
      arr[i] = cur->p;
   }
   assert(i == nn); // making sure that Vector_nnz() works as it should 
   bn_initialize(arr[nn].index_bn); // the -1 termination of the array:
   bn_setl(arr[nn].index_bn, -1); // arr[nn].index = -1;
   arr[nn].index = -1; // the -1 termination of the array
   arr[nn].val = 0.0;
   return arr;
}

/*--------------------------------------------------------------------*/

/* Apply f to every component of the vector and return the result as a
   new vector. Return NULL if insufficient memory available. 
   Note: DOES NOT free or modify the original vector. */

Vector Vector_map(Vector vec, struct Pair (*f)(struct Pair, void*), 
      void* extra) {
   Vector target = Vector_new_empty();
   for (Vector cur = vec; cur != NULL; cur = cur->next) {
      struct Pair p_new = f(cur->p, extra);
      // // remove after refactoring !!
      // assert (bn_lval(p_new.index_bn) == p_new.index);
      if (bn_lval(p_new.index_bn) != p_new.index) {
         bn_setl(p_new.index_bn, p_new.index); // remove after refactoring !!
      }
      if (bn_cmpl(p_new.index_bn, -1) == 0) continue; // if (p_new.index == -1)
      assert(bn_sgn(p_new.index_bn) >= 0); // assert(p_new.index >= 0);
      Vector target_new = Vector_addto_bn(target, p_new.index_bn, p_new.val);
      if (target_new == NULL) { // insufficient memory
         Vector_free(target);
         return NULL;
      }
      target = target_new;
   }
   return target;
}

/* Same as map, but where the function may map one component to multiple.
   The Pair array returned by f must be terminated by a -1 index. 
   Note: DOES NOT free or modify the original vector. */

// Temporary version of this function. Causes memory leaks, but should
// no longer do so when the marked code is removed
Vector Vector_map_mult(Vector vec, struct Pair* (*f)(struct Pair, void*), 
      void* extra) {
   Vector target = Vector_new_empty();
   for (Vector cur = vec; cur != NULL; cur = cur->next) {
      struct Pair* p_arr = f(cur->p, extra);
      if (p_arr == NULL || p_arr->index == -1) //remove post refactoring !!
            continue; // nothing to do         //remove post refactoring !!
      if (p_arr == NULL || bn_cmpl(p_arr->index_bn, -1) == 0) 
         continue; // nothing to do
      struct Pair* p_new = p_arr;
      //remove post refactoring!!
      assert (bn_lval(p_new->index_bn) == p_new->index);
      // fprintf(stderr, "1 marker    | %ld\n", 2131L);
      while ((p_new->index != -1)         // remove 1st part post refactoring !!
            && bn_cmpl(p_new->index_bn, -1) != 0) { 
         if (bn_sgn(p_new->index_bn) == 0 //remove post refactoring !!
               || bn_cmpl(p_new->index_bn, p_new->index) != 0) {
            bn_initialize(p_new->index_bn); //remove post refactoring !!
            bn_setl(p_new->index_bn, p_new->index); //remove post refactoring !!
         }
         // bool array_misformatted = p_new->index < 0;
         bool array_misformatted = bn_sgn(p_new->index_bn) < 0;
         assert(!array_misformatted); // can be violated if array not terminated
         Vector target_new = 
               Vector_addto_bn(target, p_new->index_bn, p_new->val);
         if (target_new == NULL) { // insufficient memory
            Vector_free(target);
            return NULL;
         }
         target = target_new;
         p_new++;
      }
      // free(p_arr);
   }
   return target;
}

/*--------------------------------------------------------------------*/

// /* Deprecated: Return a reduced version of the vector with all 
//    zero and almost-zero values removed. */

// Vector Vector_reduce(Vector vec);

/*--------------------------------------------------------------------*/

/* Return a vector with all values multiplied by the scalar c.
   The _copy version preserves the argument vector so that 
   both the argument and result vector are in memory. */
   
Vector Vector_scalar_prod_copy(const Vector vec, double complex c) {
   struct Pair prod(struct Pair p, void* extra) {
      p.val = c * p.val;
      return p;
   }
   return Vector_map(vec,&prod,NULL);
}

/* Return the vector with all values multiplied by the scalar c. */

Vector Vector_scalar_prod(Vector vec, double complex c) {
   Vector target = Vector_scalar_prod_copy(vec,c);
   Vector_free(vec);
   return target;
}

/*--------------------------------------------------------------------*/

/* Return the vector with values that are the complex conjugates of 
   those of vec. 
   The _copy version preserves the argument vector so that 
   both the argument and result vector are in memory. */

Vector Vector_conj_copy(const Vector vec) {
   struct Pair conjug(struct Pair p, void* extra) {
      p.val = conj(p.val);
      return p;
   }
   return Vector_map(vec,&conjug,NULL);
}

/* Return the vector with values that are the complex conjugates of 
   those of vec. */

Vector Vector_conj(Vector vec) {
   Vector target = Vector_conj_copy(vec);
   Vector_free(vec);
   return target;
}

/*--------------------------------------------------------------------*/

/* Return a vector that is the sum of the two vectors. 
   Does not modify either argument vector. */

Vector Vector_sum(const Vector v1, const Vector v2) {
   Vector target = Vector_new_empty();
   for (Vector cur = v1; cur != NULL; cur = cur->next) {
      Vector target_new = Vector_addto_bn(target, cur->p.index_bn, cur->p.val);
      if (target_new == NULL) { // insufficient memory
         Vector_free(target);
         return NULL;
      }
      target = target_new;
   }
   for (Vector cur = v2; cur != NULL; cur = cur->next) {
      Vector target_new = Vector_addto_bn(target, cur->p.index_bn, cur->p.val);
      if (target_new == NULL) { // insufficient memory
         Vector_free(target);
         return NULL;
      }
      target = target_new;
   }
   return target;
}

/*--------------------------------------------------------------------*/

/* Return a vector that is the element-wise product of the two vectors. 
   Does not modify either argument vector. */

Vector Vector_element_prod(const Vector v1, const Vector v2) {
   Vector target = Vector_new_empty();
   Vector cur1 = v1;
   Vector cur2 = v2;
   while (cur1 && cur2) {
      // while (cur1 != NULL && cur1->p.index < cur2->p.index) {
      while (cur1 != NULL && bn_cmp(cur1->p.index_bn, cur2->p.index_bn) < 0) {
         cur1 = cur1->next;
      }
      if (cur1 == NULL) break;
      // while (cur2 != NULL && cur2->p.index < cur1->p.index) {
      while (cur2 != NULL && bn_cmp(cur2->p.index_bn, cur1->p.index_bn) < 0) {
         cur2 = cur2->next;
      }
      if (cur2 == NULL) break;
      // if (cur2->p.index == cur1->p.index) {
      if (bn_cmp(cur2->p.index_bn, cur1->p.index_bn) == 0) {
         double complex new_val = cur1->p.val * cur2->p.val;
         Vector target_new = 
               Vector_insert_bn(target, cur1->p.index_bn, new_val);
         if (target_new == NULL) { // insufficient memory
            Vector_free(target);
            return NULL;
         }
         target = target_new;
         cur1 = cur1->next;
         cur2 = cur2->next;
      }
   }
   return target;
}

/*--------------------------------------------------------------------*/

/* Return the inner product of the two vectors. 
   Does not modify either argument vector. */

double complex Vector_inner_prod(const Vector v1, const Vector v2) {
   Vector v1_c = Vector_conj_copy(v1);
   Vector x = Vector_element_prod(v1_c,v2);
   double complex total = 0.0;
   for (Vector cur = x; cur != NULL; cur = cur->next) {
      total += cur->p.val;
   }
   Vector_free(v1_c);
   return total;
}

/*--------------------------------------------------------------------*/

/* Return the L2 norm or magnitude of a vector.
   Does not modify the argument vector. */

double complex Vector_norm(const Vector vec) {
   double complex total = 0.0;
   for (Vector cur = vec; cur != NULL; cur = cur->next) {
      total += conj(cur->p.val) * cur->p.val;
   }
   return sqrt(total);
   // // alternatively
   // return sqrt(Vector_inner_prod(vec, vec));
}



/*--------------------------------------------------------------------*/
/* Unit testing for vector.c                                          */
/*--------------------------------------------------------------------*/

bool isValid(const Vector vec) {
   // if vec is empty, then it is valid
   if (vec == NULL) return true;
   
   // check for no cycles
   Vector tort = vec;
   Vector hare = vec->next; // vec != NULL because of check above
   while (tort != NULL && hare != NULL) {
      if (tort == hare) // found a cycle
         return false;
      tort = tort->next; // tortoise moves one forward
      if (hare->next == NULL) break; // must have reached the end
      hare = hare->next->next; // hare moves two foreward
   }
   
   // check for strictly ascending order. Also detects duplicates
   Vector cur = vec;
   Vector nex = vec->next; // vec != NULL because of check above
   while (nex != NULL) {
      // if (nex->p.index < 0 || cur->p.index < 0) // negative indices
      if (bn_sgn(nex->p.index_bn) < 0 || bn_sgn(cur->p.index_bn) < 0) 
         return false;                    // negative indices
      // if (nex->p.index <= cur->p.index) // indices in wrong order
      if (bn_cmp(nex->p.index_bn, cur->p.index_bn) <= 0) 
         return false;                    // indices in wrong order
      cur = nex;
      nex = cur->next;
   }
   
   // Remove after refactor !!! uses wrong type of index
   // check for strictly ascending order. Also detects duplicates
   cur = vec;
   nex = vec->next; // vec != NULL because of check above
   while (nex != NULL) {
      if (nex->p.index < 0 || cur->p.index < 0) // negative indices
         return false;
      if (nex->p.index <= cur->p.index) // indices in wrong order
         return false;
      cur = nex;
      nex = cur->next;
   }
   
   // everything checks out!
   return true;
}

/* Prints out the Pairs contained in the vector, one per line. */

void Vector_print(const Vector vec) {
   struct Pair p_end;
   bn_initialize(p_end.index_bn);
   bn_setl(p_end.index_bn, -1);
   struct Pair print_pair (struct Pair p, void* extra) {
      bn_print(p.index_bn);
      printf(": %.1f+%.1fi\n", creal(p.val), cimag(p.val));
      return p_end;
   }
   printf("NNZ: %ld\n", Vector_nnz(vec));
   Vector_free(Vector_map(vec,&print_pair,NULL)); 
   fflush(stdout);
   bn_clear(p_end.index_bn);
}

/* Functions that return whether or not a Vector and an array of 
   Pairs contain the same set of elements. 
   The argument len should specify the length of arr. */
bool vector_contains_array(const Vector v, struct Pair* arr, long len) {
   for (long i = 0; i < len; i++) {
      if (Vector_get(v,arr[i].index) != arr[i].val)
         return false;
   }
   return true;
}
bool vector_contained_in_array(Vector v, struct Pair* arr, long len) {
   // make a copy of arr
   struct Pair a [len+1];
   for (long i = 0; i < len; i++) {
      a[i] = arr[i];
   }
   a[len].index = -1; // -1 indicates okay, -2 indicates a problem
   struct Pair* check_in_arr (struct Pair p, void* extra) {
      bool found = false;
      for (long i = 0; i < len; i++) {
         if (a[i].index == p.index) {
            if (a[i].val != p.val) 
               a[len].index = -2;
            found = true;
            break;
         }
      }
      if (!found) a[len].index = -2;
      return NULL;
   }
   Vector nothing = Vector_map_mult(v,&check_in_arr,NULL);
   return (a[len].index == -1);
}
bool vector_equals_array(Vector v, struct Pair* arr, long len) {
   return vector_contains_array(v,arr,len) 
         && vector_contained_in_array(v,arr,len);
}

void test_constructors (int n) {
   assert(n>=4);
   Vector v;
   struct Pair p; // arbitrary sample pair
   p.index = 382;
   p.val = 2.4 - 7.3 * I;
   // test the empty constructor
   v = Vector_new_empty();
   assert(isValid(v));
   assert(Vector_isEmpty(v));
   assert(Vector_nnz(v) == 0);
   assert(Vector_get(v,0) == 0.0);
   Vector_free(v);
   // test the "singleton" constructor
   v = Vector_new(p.index, p.val);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == 1);
   assert(Vector_get(v,p.index+1) == 0.0);
   assert(Vector_get(v,p.index) == p.val);
   Vector_free(v);
}
void test_insert (int n) {
   assert(n>=4);
   Vector v;
   struct Pair a[n+1];
   struct Pair p; // arbitrary sample pair
   p.index = 382 + n;
   p.val = 2.4 - 7.3 * I;
   
   // insert into an empty vector
   v = Vector_new_empty();
   assert(isValid(v));
   assert(Vector_isEmpty(v));
   assert(Vector_nnz(v) == 0);
   assert(Vector_get(v,p.index) == 0.0);
   v = Vector_insert(v, p.index, p.val);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == 1);
   assert(Vector_get(v,p.index) == p.val);
   Vector_free(v);
   // insert into a small vector already containing the index
   v = Vector_new(p.index,1.0);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == 1);
   assert(Vector_get(v,p.index) == 1.0);
   v = Vector_insert(v, p.index, p.val);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == 1);
   assert(Vector_get(v,p.index) == p.val);
   Vector_free(v);
   // insert into a small vector not containing the index
   v = Vector_new(p.index+1,1.0);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == 1);
   assert(Vector_get(v,p.index) == 0.0);
   v = Vector_insert(v, p.index, p.val);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == 2);
   assert(Vector_get(v,p.index) == p.val);
   Vector_free(v);
   
   // insert many values in ascending order
   v = Vector_new_empty();
   for (int i = 0; i < n; i++) {
      a[i].index = i;
      a[i].val = i*i - i * I + 1;
      v = Vector_insert(v,a[i].index,a[i].val);
      assert(v); // check that did not run out of memory
   }
   assert(isValid(v));
   assert(vector_equals_array(v,a,n));
   assert(Vector_nnz(v) == n);
   Vector_free(v);
   
   // insert in descending order
   v = Vector_new_empty();
   for (int i = n-1; i >= 0; i--) {
      a[i].index = i;
      a[i].val = i*i - i * I + 1;
      v = Vector_insert(v,a[i].index,a[i].val);
      assert(v);
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(vector_equals_array(v,a,n));
   assert(Vector_nnz(v) == n);
   Vector_free(v);
   
   // insert in random order
   v = Vector_new_empty();
   for (int i = n-1; i >= 0; i--) {
      do {
         a[i].index = rand() % (n*10);
      } while (Vector_get(v,a[i].index) != 0.0);
      a[i].val = i*i - i * I + 1;
      v = Vector_insert(v,a[i].index,a[i].val);
      assert(v);
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(vector_equals_array(v,a,n));
   assert(Vector_nnz(v) == n);
   Vector_free(v);
}
void test_addto (int n) {
   assert(n>=4);
   Vector v;
   struct Pair a[n+1];
   struct Pair p; // arbitrary sample pair
   p.index = 382 + n;
   p.val = 2.4 - 7.3 * I;
   long new_index;
   
   // add to an empty vector
   v = Vector_new_empty();
   assert(isValid(v));
   assert(Vector_isEmpty(v));
   assert(Vector_nnz(v) == 0);
   assert(Vector_get(v,p.index) == 0.0);
   v = Vector_addto(v, p.index, p.val);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == 1);
   assert(Vector_get(v,p.index) == p.val);
   Vector_free(v);
   // add to a small vector already containing the index
   v = Vector_new(p.index,1.0);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == 1);
   assert(Vector_get(v,p.index) == 1.0);
   v = Vector_addto(v, p.index, p.val);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == 1);
   assert(Vector_get(v,p.index) == 1.0 + p.val);
   Vector_free(v);
   // add to a small vector not containing the index
   v = Vector_new(p.index+1,1.0);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == 1);
   assert(Vector_get(v,p.index) == 0.0);
   v = Vector_addto(v, p.index, p.val);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == 2);
   assert(Vector_get(v,p.index) == p.val);
   Vector_free(v);
   
   // add to a large vector already containing the index
   v = Vector_new_empty();
   for (int i = n-1; i >= 0; i--) {
      do {
         a[i].index = rand() % (n*10);
      } while (Vector_get(v,a[i].index) != 0.0);
      a[i].val = i*i - i * I + 1;
      v = Vector_insert(v,a[i].index,a[i].val);
      assert(v);
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(vector_equals_array(v,a,n));
   assert(Vector_nnz(v) == n);
   new_index = a[0].index;
   assert(Vector_get(v,new_index) == a[0].val);
   v = Vector_addto(v, new_index, p.val); // <<<<
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == n);
   assert(Vector_get(v,new_index) == a[0].val + p.val);
   Vector_free(v);
   // add to a large vector not containing the index
   v = Vector_new_empty();
   for (int i = n-1; i >= 0; i--) {
      do {
         a[i].index = rand() % (n*10);
      } while (Vector_get(v,a[i].index) != 0.0);
      a[i].val = i*i - i * I + 1;
      v = Vector_insert(v,a[i].index,a[i].val);
      assert(v);
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(vector_equals_array(v,a,n));
   assert(Vector_nnz(v) == n);
   do { // let new_index be new random index not in the vector
      new_index = rand() % (n*10);
   } while (Vector_get(v,new_index) != 0.0);
   assert(Vector_get(v,new_index) == 0.0);
   v = Vector_addto(v, new_index, p.val); // <<<<
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == n+1);
   assert(Vector_get(v,new_index) == p.val);
   Vector_free(v);
   
   // add many values in ascending order
   v = Vector_new_empty();
   for (int i = 0; i < n; i++) {
      a[i].index = i;
      a[i].val = i*i - i * I + 1;
      v = Vector_addto(v,a[i].index,a[i].val);
      assert(v); // check that did not run out of memory
   }
   assert(isValid(v));
   assert(vector_equals_array(v,a,n));
   assert(Vector_nnz(v) == n);
   Vector_free(v);
   
   // add values in descending order
   v = Vector_new_empty();
   for (int i = n-1; i >= 0; i--) {
      a[i].index = i;
      a[i].val = i*i - i * I + 1;
      v = Vector_addto(v,a[i].index,a[i].val);
      assert(v);
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(vector_equals_array(v,a,n));
   assert(Vector_nnz(v) == n);
   Vector_free(v);
   
   // add values in random order
   v = Vector_new_empty();
   for (int i = 0; i < n; i++) {
      do {
         a[i].index = rand() % (n*10);
      } while (Vector_get(v,a[i].index) != 0.0);
      a[i].val = i*i - i * I + 1;
      v = Vector_addto(v,a[i].index,a[i].val);
      assert(v);
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(vector_equals_array(v,a,n));
   assert(Vector_nnz(v) == n);
   Vector_free(v);
}
void test_get (int n) {
   assert(n>=4);
   Vector v;
   struct Pair a[n+1];
   struct Pair p; // arbitrary sample pair
   p.index = 382 + n;
   p.val = 2.4 - 7.3 * I;
   long new_index;
   
   // get from an empty vector
   v = Vector_new_empty();
   assert(isValid(v));
   assert(Vector_isEmpty(v));
   assert(Vector_nnz(v) == 0);
   for (int i = 0; i < n; i++) {
      assert(Vector_get(v,i) == 0.0);
      // v should not be changed
      assert(isValid(v));
      assert(Vector_isEmpty(v));
      assert(Vector_nnz(v) == 0);
   }
   Vector_free(v);
   
   // get from a singleton vector
   long ind = p.index+n+2;
   v = Vector_new(ind,p.val);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == 1);
   for (int i = 0; i < n; i++) {
      assert(Vector_get(v,i) == 0.0);
      // v should not be changed
      assert(isValid(v));
      assert(!Vector_isEmpty(v));
      assert(Vector_nnz(v) == 1);
   }
   assert(Vector_get(v,ind) == p.val);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == 1);
   Vector_free(v);
   
   // get from a large vector
   v = Vector_new_empty();
   for (int i = 0; i < n; i++) {
      a[i].index = i+n;
      a[i].val = i*i - i * I + 1;
      v = Vector_addto(v,a[i].index,a[i].val);
      assert(v); // check that did not run out of memory
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(vector_equals_array(v,a,n));
   assert(Vector_nnz(v) == n);
   for (int i = 0; i < n; i++) { // indices not in the vector
      assert(Vector_get(v,i) == 0.0);
      // v should not be changed
      assert(isValid(v));
      assert(!Vector_isEmpty(v));
      assert(vector_equals_array(v,a,n));
      assert(Vector_nnz(v) == n);
   }
   for (int i = 0; i < n; i++) { // indices in the vector
      assert(Vector_get(v,a[i].index) == a[i].val);
      // v should not be changed
      assert(isValid(v));
      assert(!Vector_isEmpty(v));
      assert(vector_equals_array(v,a,n));
      assert(Vector_nnz(v) == n);
   }
}
void test_delete (int n) {
   assert(n>=4);
   Vector v;
   struct Pair a[n+1];
   struct Pair p; // arbitrary sample pair
   p.index = 382 + n;
   p.val = 2.4 - 7.3 * I;
   
   // delete from an empty vector
   v = Vector_new_empty();
   assert(isValid(v));
   assert(Vector_isEmpty(v));
   assert(Vector_nnz(v) == 0);
   assert(Vector_get(v,p.index) == 0.0);
   v = Vector_delete(v, p.index);
   assert(isValid(v));
   assert(Vector_isEmpty(v));
   assert(Vector_nnz(v) == 0);
   assert(Vector_get(v,p.index) == 0.0);
   Vector_free(v);
   // delete from a singleton vector
   v = Vector_new(p.index,p.val);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == 1);
   assert(Vector_get(v,p.index) == p.val);
   v = Vector_delete(v, p.index);
   assert(isValid(v));
   assert(Vector_isEmpty(v)); // should now be empty
   assert(Vector_nnz(v) == 0);
   assert(Vector_get(v,p.index) == 0.0);
   Vector_free(v);
   // delete from a large vector containing the index
   v = Vector_new_empty();
   for (int i = 0; i < n; i++) {
      do {
         a[i].index = rand() % (n*10);
      } while (Vector_get(v,a[i].index) != 0.0);
      a[i].val = i*i - i * I + 1;
      v = Vector_insert(v,a[i].index,a[i].val);
      assert(v);
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(vector_equals_array(v,a,n));
   assert(Vector_nnz(v) == n);
   long ind = a[n-1].index; // select the last (nth) index from a
   assert(Vector_get(v,ind) == a[n-1].val);
   v = Vector_delete(v, ind); // <<<<
   assert(isValid(v));
   assert(Vector_get(v,ind) == 0.0); // it should now be gone
   assert(!Vector_isEmpty(v)); // as long as n > 1
   assert(Vector_nnz(v) == n - 1);
   assert(vector_equals_array(v,a,n-1)); // equivalent the first n-1 pairs of a
   Vector_free(v);
   // delete from a large vector not containing the index
   v = Vector_new_empty();
   for (int i = 0; i < n; i++) {
      do {
         a[i].index = rand() % (n*10);
      } while (Vector_get(v,a[i].index) != 0.0);
      a[i].val = i*i - i * I + 1;
      v = Vector_insert(v,a[i].index,a[i].val);
      assert(v);
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(vector_equals_array(v,a,n));
   assert(Vector_nnz(v) == n);
   do { // select a new random index not in the vector
      ind = rand() % (n*10);
   } while (Vector_get(v,ind) != 0.0);
   assert(Vector_get(v,ind) == 0.0);
   v = Vector_delete(v, ind); // <<<< This should have no effect
   assert(isValid(v));
   assert(Vector_get(v,ind) == 0.0);
   assert(!Vector_isEmpty(v));
   assert(vector_equals_array(v,a,n));
   assert(Vector_nnz(v) == n);
   Vector_free(v);
   // delete all indices from a vector
   v = Vector_new_empty();
   for (int i = 0; i < n; i++) {
      do {
         a[i].index = rand() % (n*10);
      } while (Vector_get(v,a[i].index) != 0.0);
      a[i].val = i*i - i * I + 1;
      v = Vector_insert(v,a[i].index,a[i].val);
      assert(v);
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(vector_equals_array(v,a,n));
   assert(Vector_nnz(v) == n);
   for (int i = 0; i < n; i++) {
      assert(Vector_get(v,a[i].index) == a[i].val);
   }
   for (int i = 0; i < n; i++) { // <<<<
      assert(isValid(v));
      assert(Vector_nnz(v) == n-i);
      v = Vector_delete(v,a[i].index);
   }
   assert(isValid(v));
   assert(Vector_isEmpty(v)); // v should now be empty
   assert(Vector_nnz(v) == 0);
   for (int i = 0; i < n; i++) {
      assert(Vector_get(v,a[i].index) == 0.0);
   }
   Vector_free(v);
}
void test_empty (int n) {
   assert(n>=4);
   Vector v;
   struct Pair p;
   p.index = 382;
   p.val = 2.4 - 7.3 * I;
   // the empty vector should be empty
   v = Vector_new_empty();
   assert(Vector_isEmpty(v));
   // but inserting an element should make it nonempty
   v = Vector_insert(v, p.index, p.val);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   Vector_free(v);
   // the "singleton" vector should be nonempty
   v = Vector_new(p.index, p.val);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   // but deleting its lone element should make it empty
   v = Vector_delete(v, p.index);
   assert(isValid(v));
   assert(Vector_isEmpty(v));
   Vector_free(v);
}
void test_nnz (int n) {
   // for brevity, the comments here refer to the number 
   // of nonzeroes as the "length" of the vector
   assert(n>=4);
   Vector v;
   struct Pair a[n+1];
   struct Pair p; // arbitrary sample pair
   p.index = 382 + n;
   p.val = 2.4 - 7.3 * I;
   
   // the empty vector should have length 0
   v = Vector_new_empty();
   assert(Vector_nnz(v) == 0);
   // but inserting an element should give it length 1
   v = Vector_insert(v, p.index, p.val);
   assert(isValid(v));
   assert(Vector_nnz(v) == 1);
   Vector_free(v);
   // the "singleton" vector should have length 1
   v = Vector_new(p.index, p.val);
   assert(isValid(v));
   assert(Vector_nnz(v) == 1);
   // but deleting its lone element should give it length 0
   v = Vector_delete(v, p.index);
   assert(isValid(v));
   assert(Vector_nnz(v) == 0);
   Vector_free(v);
   // a vector of n numbers should have length n
   v = Vector_new_empty();
   for (int i = 0; i < n; i++) {
      do {
         a[i].index = rand() % (n*10);
      } while (Vector_get(v,a[i].index) != 0.0);
      a[i].val = i*i - i * I + 1;
      v = Vector_insert(v,a[i].index,a[i].val);
      assert(v);
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(vector_equals_array(v,a,n));
   assert(Vector_nnz(v) == n); // <<<<
   // deleting two elements should give it length n - 2
   long ind1 = a[n-1].index; // select the last two indices from a
   long ind2 = a[n-2].index;
   v = Vector_delete(Vector_delete(v, ind2), ind1);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(vector_equals_array(v,a,n-2));
   assert(Vector_nnz(v) == n-2); // <<<<
   // but inserting an element should give it length n - 1
   v = Vector_insert(v, ind2, a[n-2].val);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(vector_equals_array(v,a,n-1));
   assert(Vector_nnz(v) == n-1); // <<<<
   // adding to an existing index should not change the length
   v = Vector_addto(v, ind2, 1.0);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == n-1); // <<<<
   Vector_free(v);
}
void test_contents (int n) {
   assert(n>=4);
   Vector v;
   struct Pair a[n+1];
   struct Pair p; // arbitrary sample pair
   p.index = 382 + n;
   p.val = 2.4 - 7.3 * I;
   struct Pair* arr;
   
   // contents of an empty vector
   v = Vector_new_empty();
   assert(isValid(v));
   assert(Vector_isEmpty(v));
   assert(Vector_nnz(v) == 0);
   arr = Vector_contents(v);
   assert(arr[0].index == -1); // empty but terminated by a -1
   free(arr);
   Vector_free(v);
   // contents of a singleton vector
   v = Vector_new(p.index, p.val);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == 1);
   arr = Vector_contents(v); // <<<<
   assert(arr[0].index == p.index);
   assert(arr[0].val == p.val);
   assert(arr[1].index == -1);
   assert(vector_equals_array(v,arr,1));
   free(arr);
   Vector_free(v);
   // contents of a large vector
   v = Vector_new_empty();
   for (int i = 0; i < n; i++) {
      do {
         a[i].index = rand() % (n*10);
      } while (Vector_get(v,a[i].index) != 0.0);
      a[i].val = i*i - i * I + 1;
      v = Vector_insert(v,a[i].index,a[i].val);
      assert(v);
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(vector_equals_array(v,a,n));
   assert(Vector_nnz(v) == n);
   arr = Vector_contents(v); // <<<<
   assert(arr[n].index == -1);
   assert(vector_equals_array(v,arr,n));
   free(arr);
   Vector_free(v);
}
void test_map (int n) {
   assert(n>=4);
   Vector v;
   Vector v_old;
   struct Pair a[n+1];
   struct Pair b[2*n+1];
   struct Pair p; // arbitrary sample pair
   p.index = 382 + n;
   p.val = 2.4 - 7.3 * I;
   
   // use map to negate all the values
   v = Vector_new_empty();
   for (int i = 0; i < n; i++) {
      do {
         a[i].index = rand() % (n*10);
      } while (Vector_get(v,a[i].index) != 0.0);
      a[i].val = i*i - i * I + 1;
      b[i].index = a[i].index;
      b[i].val = -1 * a[i].val;
      v = Vector_insert(v,a[i].index,a[i].val);
      assert(v);
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == n);
   assert(vector_equals_array(v,a,n)); // initially like a
   struct Pair negate(struct Pair p, void* extra) {
      p.val = -1 * p.val;
      return p;
   }
   v_old = v;
   v = Vector_map(v,&negate,NULL); // <<<
   Vector_free(v_old);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == n);
   assert(vector_equals_array(v,b,n)); // now negated like b
   Vector_free(v);
   
   // use map to flip the indices i -> 10 * n - i
   v = Vector_new_empty();
   for (int i = 0; i < n; i++) {
      do {
         a[i].index = rand() % (n*10);
      } while (Vector_get(v,a[i].index) != 0.0);
      a[i].val = i*i - i * I + 1;
      b[i].index = 10 * n - a[i].index;
      b[i].val = a[i].val;
      v = Vector_insert(v,a[i].index,a[i].val);
      assert(v);
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == n);
   assert(vector_equals_array(v,a,n)); // initially like a
   struct Pair p_new;
   bn_initialize(p_new.index_bn);
   struct Pair flip(struct Pair p, void* extra) {
      bn_set(p_new.index_bn, p.index_bn);
      bn_sub_from(p_new.index_bn, 10 * n);
      bn_negate(p_new.index_bn);
      p_new.index = 10 * n - p.index;
      p_new.val = p.val;
      return p_new;
   }
   v_old = v;
   v = Vector_map(v,&flip,NULL); // <<<
   bn_clear(p_new.index_bn);
   Vector_free(v_old);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == n);
   assert(vector_equals_array(v,b,n)); // now flipped like b
   Vector_free(v);
   
   // map each index to itself and the one after it
   v = Vector_new_empty();
   for (int i = 0; i < n; i++) {
      do {
         a[i].index = rand() % (n*20) + 1;
      } while ((Vector_get(v,a[i].index) != 0.0) // make sure indices are
            || (Vector_get(v,a[i].index + 1) != 0.0) // not adjacent
            || (Vector_get(v,a[i].index - 1) != 0.0));
      a[i].val = i*i - i * I + 1;
      b[2*i].index = a[i].index;
      b[2*i].val = a[i].val;
      b[2*i+1].index = a[i].index + 1;
      b[2*i+1].val = a[i].val;
      v = Vector_insert(v,a[i].index,a[i].val);
      assert(v);
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == n);
   assert(vector_equals_array(v,a,n)); // initially like a
   struct Pair d[3]; // the array to be returned by duplic
   bn_initialize(d[0].index_bn);
   bn_initialize(d[1].index_bn);
   bn_initialize(d[2].index_bn);
   struct Pair* duplic(struct Pair p, void* extra) {
      bn_setl(d[0].index_bn, p.index);
      d[0].index = p.index;
      d[0].val = p.val;
      bn_setl(d[1].index_bn, p.index + 1);
      d[1].index = p.index + 1;
      d[1].val = p.val;
      bn_setl(d[2].index_bn, -1);
      d[2].index = -1;
      return d;
   }
   v_old = v;
   v = Vector_map_mult(v,&duplic,NULL); // <<<
   bn_clear(d[0].index_bn);
   bn_clear(d[1].index_bn);
   bn_clear(d[2].index_bn);
   Vector_free(v_old);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == 2 * n);
   assert(vector_equals_array(v,b,2*n)); // now duplicated like b
   Vector_free(v);
   
   // map each single index to all indices
   v = Vector_new_empty();
   a[0].index = 0;
   a[0].val = 1.0;
   a[1].index = 10;
   a[1].val = 2.0;
   a[2].index = 53;
   a[2].val = 3.0;
   for (int i = 0; i < 3; i++) {
      v = Vector_insert(v,a[i].index,a[i].val);
      assert(v);
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == 3);
   assert(vector_equals_array(v,a,3)); // initially like a
   struct Pair e[n+1]; // the array to be returned by to_all
   for (int i = 0; i < n+1; i++) {
      bn_initialize(e[i].index_bn);
   }
   struct Pair* to_all(struct Pair p, void* extra) {
      for (int i = 0; i < n; i++) {
         bn_setl(e[i].index_bn, i);
         e[i].index = i;
         e[i].val = p.val;
      }
      e[n].index = -1;
      return e;
   }
   v_old = v;
   v = Vector_map_mult(v,&to_all,NULL); // <<<
   for (int i = 0; i < n+1; i++) {
      bn_clear(e[i].index_bn);
   }
   Vector_free(v_old);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == n);
   for (int i = 0; i < n; i++) {
      b[i].index = i;
      b[i].val = 1.0 + 2.0 + 3.0;
   }
   assert(vector_equals_array(v,b,n)); // now evenly distributed like b
   Vector_free(v);
}
void test_scalar_prod (int n) {
   assert(n>=4);
   Vector v;
   struct Pair a[n+1];
   struct Pair b[n+1];
   struct Pair p; // arbitrary sample pair
   p.index = 382 + n;
   p.val = 2.4 - 7.3 * I;
   
   // test multiplying by I
   v = Vector_new_empty();
   for (int i = 0; i < n; i++) {
      do {
         a[i].index = rand() % (n*10);
      } while (Vector_get(v,a[i].index) != 0.0);
      a[i].val = i*i - i * I + 1;
      b[i].index = a[i].index;
      b[i].val = I * a[i].val;
      v = Vector_insert(v,a[i].index,a[i].val);
      assert(v);
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == n);
   assert(vector_equals_array(v,a,n));
   v = Vector_scalar_prod(v,I); // <<<< multiply by I = sqrt(-1)
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == n);
   assert(vector_equals_array(v,b,n));
   v = Vector_scalar_prod(v,I); // <<<< multiply by I three more times
   v = Vector_scalar_prod(v,I);
   v = Vector_scalar_prod(v,I);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == n);
   assert(vector_equals_array(v,a,n)); // should bring it back to what it was
   Vector_free(v);
   // test multiplying by 0
   v = Vector_new_empty();
   for (int i = 0; i < n; i++) {
      do {
         a[i].index = rand() % (n*10);
      } while (Vector_get(v,a[i].index) != 0.0);
      a[i].val = i*i - i * I + 1;
      b[i].index = a[i].index;
      b[i].val = 0.0;
      v = Vector_insert(v,a[i].index,a[i].val);
      assert(v);
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == n);
   assert(vector_equals_array(v,a,n));
   v = Vector_scalar_prod(v,0.0); // <<<<
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == n);
   assert(vector_equals_array(v,b,n));
      // assert(Vector_nnz(Vector_reduce(v)) == 0);
   Vector_free(v);
   // test dividing by two
   v = Vector_new_empty();
   for (int i = 0; i < n; i++) {
      do {
         a[i].index = rand() % (n*10);
      } while (Vector_get(v,a[i].index) != 0.0);
      a[i].val = i*i - i * I + 1;
      b[i].index = a[i].index;
      b[i].val = 0.5 * a[i].val;
      v = Vector_insert(v,a[i].index,a[i].val);
      assert(v);
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == n);
   assert(vector_equals_array(v,a,n));
   v = Vector_scalar_prod(v,0.5); // <<<<
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == n);
   assert(vector_equals_array(v,b,n));
   Vector_free(v);
}
void test_sum (int n) {
   assert(n>=4);
   Vector t;
   Vector u;
   Vector v;
   struct Pair a[n+1];
   struct Pair b[n+1];
   struct Pair c[2*n+1];
   struct Pair p; // arbitrary sample pair
   p.index = 382 + n;
   p.val = 2.4 - 7.3 * I;
   
   // sum two empty vectors
   t = Vector_new_empty();
   u = Vector_new_empty();
   assert(isValid(t));
   assert(Vector_isEmpty(t));
   assert(Vector_nnz(t) == 0);
   assert(isValid(u));
   assert(Vector_isEmpty(u));
   assert(Vector_nnz(u) == 0);
   v = Vector_sum(t,u); // <<<<
   assert(isValid(v));
   assert(Vector_isEmpty(v));
   assert(Vector_nnz(v) == 0);
   Vector_free(t);
   Vector_free(u);
   Vector_free(v);
   // sum an empty vector and a nonempty vector
   u = Vector_new_empty();
   v = Vector_new_empty();
   for (int i = 0; i < n; i++) {
      do {
         a[i].index = rand() % (n*10);
      } while (Vector_get(v,a[i].index) != 0.0);
      a[i].val = i*i - i * I + 1;
      b[i].index = a[i].index;
      b[i].val = 0.5 * a[i].val;
      v = Vector_insert(v,a[i].index,a[i].val);
      assert(v);
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == n);
   assert(vector_equals_array(v,a,n));
   t = Vector_sum(u,v); // <<<< t should be equal to v since u == 0
   assert(isValid(t));
   assert(!Vector_isEmpty(t));
   assert(Vector_nnz(t) == n);
   assert(vector_equals_array(t,a,n));
   Vector_free(t);
   Vector_free(u);
   Vector_free(v);
   // sum a vector and itself
   v = Vector_new_empty();
   for (int i = 0; i < n; i++) {
      do {
         a[i].index = rand() % (n*10);
      } while (Vector_get(v,a[i].index) != 0.0);
      a[i].val = i*i - i * I + 1;
      b[i].index = a[i].index;
      b[i].val = 2 * a[i].val;
      v = Vector_insert(v,a[i].index,a[i].val);
      assert(v);
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == n);
   assert(vector_equals_array(v,a,n));
   t = v;
   v = Vector_sum(v,v); // <<<<
   Vector_free(t);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == n);
   assert(vector_equals_array(v,b,n));
   Vector_free(v);
   // sum two randomly generated vectors
   t = Vector_new_empty();
   u = Vector_new_empty();
   long ia = -1;
   long ib = -1;
   long ic = -1;
   for (int i = 0; i < n; i++) {
      if (i % 2 && i % 3) {
         a[++ia].index = i;
         b[++ib].index = i;
         c[++ic].index = i;
         a[ia].val = rand() % 50 + rand() % 50 * I + 1;
         b[ib].val = rand() % 50 + rand() % 50 * I + 1;
         c[ic].val = a[ia].val + b[ib].val;
         t = Vector_insert(t,a[ia].index,a[ia].val);
         u = Vector_insert(u,b[ib].index,b[ib].val);
         assert(t && u);
      }
      else if (i % 2) {
         a[++ia].index = i;
         c[++ic].index = i;
         a[ia].val = rand() % 50 + rand() % 50 * I + 1;
         c[ic].val = a[ia].val;
         t = Vector_insert(t,a[ia].index,a[ia].val);
         assert(t);
      }
      else if (i % 3) {
         b[++ib].index = i;
         c[++ic].index = i;
         b[ib].val = rand() % 50 + rand() % 50 * I + 1;
         c[ic].val = b[ib].val;
         u = Vector_insert(u,b[ib].index,b[ib].val);
         assert(u);
      }
   }
   assert(isValid(t));
   assert(!Vector_isEmpty(t));
   assert(Vector_nnz(t) == ia+1);
   assert(vector_equals_array(t,a,ia+1));
   assert(isValid(u));
   assert(!Vector_isEmpty(u));
   assert(Vector_nnz(u) == ib+1);
   assert(vector_equals_array(u,b,ib+1));
   v = Vector_sum(t,u); // <<<<
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == ic+1);
   assert(vector_equals_array(v,c,ic+1));
   Vector_free(t);
   Vector_free(u);
   Vector_free(v);
}
void test_element_prod (int n) {
   assert(n>=4);
   Vector t;
   Vector u;
   Vector v;
   struct Pair a[n+1];
   struct Pair b[n+1];
   struct Pair c[n+1];
   struct Pair p; // arbitrary sample pair
   p.index = 382 + n;
   p.val = 2.4 - 7.3 * I;
   
   // product of two empty vectors
   t = Vector_new_empty();
   u = Vector_new_empty();
   assert(isValid(t));
   assert(Vector_isEmpty(t));
   assert(Vector_nnz(t) == 0);
   assert(isValid(u));
   assert(Vector_isEmpty(u));
   assert(Vector_nnz(u) == 0);
   v = Vector_element_prod(t,u); // <<<<
   assert(isValid(v));
   assert(Vector_isEmpty(v));
   assert(Vector_nnz(v) == 0);
   Vector_free(t);
   Vector_free(u);
   Vector_free(v);
   // product of an empty vector and a nonempty vector
   u = Vector_new_empty();
   v = Vector_new_empty();
   for (int i = 0; i < n; i++) {
      do {
         a[i].index = rand() % (n*10);
      } while (Vector_get(v,a[i].index) != 0.0);
      a[i].val = i*i - i * I + 1;
      b[i].index = a[i].index;
      b[i].val = a[i].val;
      v = Vector_insert(v,a[i].index,a[i].val);
      assert(v);
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == n);
   assert(vector_equals_array(v,a,n));
   t = Vector_element_prod(u,v); // <<<< t should be empty (t == 0) since u == 0
   assert(isValid(t));
   assert(Vector_isEmpty(t));
   assert(Vector_nnz(t) == 0);
   Vector_free(t);
   Vector_free(u);
   Vector_free(v);
   // product of a vector and itself
   v = Vector_new_empty();
   for (int i = 0; i < n; i++) {
      do {
         a[i].index = rand() % (n*10);
      } while (Vector_get(v,a[i].index) != 0.0);
      a[i].val = i*i - i * I + 1;
      b[i].index = a[i].index;
      b[i].val = a[i].val * a[i].val;
      v = Vector_insert(v,a[i].index,a[i].val);
      assert(v);
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == n);
   assert(vector_equals_array(v,a,n));
   t = v;
   v = Vector_element_prod(v,v); // <<<<
   Vector_free(t);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == n);
   assert(vector_equals_array(v,b,n));
   Vector_free(v);
   // product of two randomly generated vectors
   t = Vector_new_empty();
   u = Vector_new_empty();
   long ia = -1;
   long ib = -1;
   long ic = -1;
   for (int i = 0; i < n; i++) {
      if (i % 2 && i % 3) {
         a[++ia].index = i;
         b[++ib].index = i;
         c[++ic].index = i;
         a[ia].val = rand() % 50 + rand() % 50 * I + 1;
         b[ib].val = rand() % 50 + rand() % 50 * I + 1;
         c[ic].val = a[ia].val * b[ib].val;
         t = Vector_insert(t,a[ia].index,a[ia].val);
         u = Vector_insert(u,b[ib].index,b[ib].val);
         assert(t && u);
      }
      else if (i % 2) {
         a[++ia].index = i;
         a[ia].val = rand() % 50 + rand() % 50 * I + 1;
         t = Vector_insert(t,a[ia].index,a[ia].val);
         assert(t);
      }
      else if (i % 3) {
         b[++ib].index = i;
         b[ib].val = rand() % 50 + rand() % 50 * I + 1;
         u = Vector_insert(u,b[ib].index,b[ib].val);
         assert(u);
      }
   }
   assert(isValid(t));
   assert(!Vector_isEmpty(t));
   assert(Vector_nnz(t) == ia+1);
   assert(vector_equals_array(t,a,ia+1));
   assert(isValid(u));
   assert(!Vector_isEmpty(u));
   assert(Vector_nnz(u) == ib+1);
   assert(vector_equals_array(u,b,ib+1));
   v = Vector_element_prod(t,u); // <<<<
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == ic+1);
   assert(vector_equals_array(v,c,ic+1));
   Vector_free(t);
   Vector_free(u);
   Vector_free(v);
}
void test_inner_prod (int n) {
   assert(n>=4);
   Vector t;
   Vector u;
   Vector v;
   struct Pair a[n+1];
   struct Pair b[n+1];
   struct Pair c[n+1];
   struct Pair p; // arbitrary sample pair
   p.index = 382 + n;
   p.val = 2.4 - 7.3 * I;
   double complex x;
   
   // product of two empty vectors
   t = Vector_new_empty();
   u = Vector_new_empty();
   assert(isValid(t));
   assert(Vector_isEmpty(t));
   assert(Vector_nnz(t) == 0);
   assert(isValid(u));
   assert(Vector_isEmpty(u));
   assert(Vector_nnz(u) == 0);
   x = Vector_inner_prod(t,u); // <<<<
   assert(x == 0.0);
   Vector_free(t);
   Vector_free(u);
   // product of an empty vector and a nonempty vector
   u = Vector_new_empty();
   v = Vector_new_empty();
   for (int i = 0; i < n; i++) {
      do {
         a[i].index = rand() % (n*10);
      } while (Vector_get(v,a[i].index) != 0.0);
      a[i].val = i*i - i * I + 1;
      b[i].index = a[i].index;
      b[i].val = a[i].val;
      v = Vector_insert(v,a[i].index,a[i].val);
      assert(v);
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == n);
   assert(vector_equals_array(v,a,n));
   x = Vector_inner_prod(u,v); // <<<< t should be empty (t == 0) since u == 0
   assert(x == 0.0);
   Vector_free(u);
   Vector_free(v);
   // product of a vector and itself
   v = Vector_new_empty();
   for (int i = 0; i < n; i++) {
      do {
         a[i].index = rand() % (n*10);
      } while (Vector_get(v,a[i].index) != 0.0);
      a[i].val = 1.0;
      b[i].index = a[i].index;
      b[i].val = a[i].val * a[i].val;
      v = Vector_insert(v,a[i].index,a[i].val);
      assert(v);
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == n);
   assert(vector_equals_array(v,a,n));
   x = Vector_inner_prod(v,v); // <<<<
   assert(x == n);
   double complex y = Vector_norm(v);
   assert(x == y*y);
   Vector_free(v);
   // product of two randomly generated vectors
   t = Vector_new_empty();
   u = Vector_new_empty();
   long ia = -1;
   long ib = -1;
   long ic = -1;
   double complex total = 0.0;
   for (int i = 0; i < n; i++) {
      if (i % 2 && i % 3) {
         a[++ia].index = i;
         b[++ib].index = i;
         c[++ic].index = i;
         a[ia].val = rand() % 50 + rand() % 50 * I + 1;
         b[ib].val = rand() % 50 + rand() % 50 * I + 1;
         c[ic].val = a[ia].val * b[ib].val;
         total += c[ic].val;
         t = Vector_insert(t,a[ia].index,a[ia].val);
         u = Vector_insert(u,b[ib].index,b[ib].val);
         assert(t && u);
      }
      else if (i % 2) {
         a[++ia].index = i;
         a[ia].val = rand() % 50 + rand() % 50 * I + 1;
         t = Vector_insert(t,a[ia].index,a[ia].val);
         assert(t);
      }
      else if (i % 3) {
         b[++ib].index = i;
         b[ib].val = rand() % 50 + rand() % 50 * I + 1;
         u = Vector_insert(u,b[ib].index,b[ib].val);
         assert(u);
      }
   }
   assert(isValid(t));
   assert(!Vector_isEmpty(t));
   assert(Vector_nnz(t) == ia+1);
   assert(vector_equals_array(t,a,ia+1));
   assert(isValid(u));
   assert(!Vector_isEmpty(u));
   assert(Vector_nnz(u) == ib+1);
   assert(vector_equals_array(u,b,ib+1));
   x = Vector_inner_prod(t,u); // <<<<
   assert(x = total);
   Vector_free(t);
   Vector_free(u);
}
void test_norm (int n) {
   assert(n>=4);
   Vector t;
   Vector u;
   Vector v;
   struct Pair a[n+1];
   struct Pair b[n+1];
   struct Pair c[n+1];
   struct Pair p; // arbitrary sample pair
   p.index = 382 + n;
   p.val = 2.4 - 7.3 * I;
   double complex x;
   
   // norm of an empty vector
   v = Vector_new_empty();
   assert(isValid(v));
   assert(Vector_isEmpty(v));
   assert(Vector_nnz(v) == 0);
   x = Vector_norm(v); // <<<<
   assert(x == 0.0);
   Vector_free(v);
   // norm of a singleton vector
   v = Vector_new(p.index,p.val);
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == 1);
   assert(Vector_get(v,p.index) == p.val);
   x = Vector_norm(v); // <<<<
   assert(x == csqrt(conj(p.val) * p.val));
   Vector_free(v);
   // norm of a randomly generated large vector
   v = Vector_new_empty();
   double complex total = 0.0;
   for (int i = 0; i < n; i++) {
      do {
         a[i].index = rand() % (n*10);
      } while (Vector_get(v,a[i].index) != 0.0);
      a[i].val = (rand() % 50 + 1) / 5 + (rand() % 50 + 1) * I / 5;
      total += conj(a[i].val) * a[i].val;
      v = Vector_insert(v,a[i].index,a[i].val);
      assert(v);
   }
   assert(isValid(v));
   assert(!Vector_isEmpty(v));
   assert(Vector_nnz(v) == n);
   assert(vector_equals_array(v,a,n));
   x = Vector_norm(v); // <<<<
   assert(creal(conj(x*x - total) * (x*x - total)) < .00001);
   double complex y = Vector_inner_prod(v,v);
   assert(creal(conj(x*x - y) * (x*x - y)) < .00001);
   Vector_free(v);
}
void test_print(int n) {
   Vector v = Vector_new_empty();
   for (int i = 0; i < n; i++) {
      v = Vector_insert(v,i,i*i - i * I + 1);
      assert(v); // check that did not run out of memory
   }
   assert(isValid(v));
   assert(Vector_nnz(v) == n);
   Vector_print(v);     // --------
   Vector_free(v);
}

/* Runs the unit testing code. Except when preblems arise, running 
   this function should produce no visible effect. */
void Vector_run_tests (int n) {
   test_constructors (n);
   test_insert (n);
   test_addto (n);
   test_get (n);
   test_delete (n);
   test_empty (n);
   test_nnz (n);
   test_contents (n);
   test_map (n);
   test_scalar_prod (n);
   test_sum (n);
   test_element_prod (n);
   test_inner_prod (n);
   test_norm (n);
   // // Uncomment this next line to test Vector_print(). Prints to stdout.
   // test_print(n);
}

int Vector_main(void) {
   srand(time(NULL));
   
   unsigned int n = 100; // Maximum size of test cases. Must be at least 4.
   Vector_run_tests (n);
   
   
   // printf("8"); fflush(stdout);
   // printf("actual: %ld, %f+%fi\n", p.index, creal(p.val), cimag(p.val));
   
   /* printf("x: %f+%fi\ntotal: %f+%fi\n", 
         creal(x),cimag(x),creal(total),cimag(total)); */
   
   // assert(false);
   
   
   // fprintf(stderr, "1 marker    | %ld\n", 2131L);
   
   return 0;
}