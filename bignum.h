/* 
 * bignum.h
 * --------
 * An set of definitions to allow conversion to use of arbitrary precision 
 * integers, which in turn allows expansions to simulations of larger numbers 
 * of qubits, only limited by system memory and time constraints.
 * 
 * Author: Barak Nehoran
 * Advisor: Iasonas Petras
 */

#ifndef BIGNUM_H
#define BIGNUM_H

#include <stdio.h>
#include <gmp.h>

// a bignum is an arbitrary precision
typedef mpz_t bignum;

// initialize the bignum
#define bn_new(n)							    \
    bignum n;                                   \
    mpz_init(n)

#define bn_initialize(n)							    \
    mpz_init(n)

// initialize and set to a long
#define bn_new_with(n, i)						\
    bignum n;                                   \
    mpz_init_set_si(n, i)

// #define bn_initialize_with(n, i)							    
//     mpz_init_set_si(n, i)

// set to a long
#define bn_setl(n, i)							\
    mpz_set_si(n, i)

// set to another bignum
#define bn_set(n, m)							\
    mpz_set(n, m)

// free the memory used by the bignum
#define bn_clear(n)							    \
    mpz_clear(n)

// return the value as a long (drops more significant bits)
#define bn_lval(n)                               \
    mpz_get_si(n)

// return the value as a double
#define bn_dval(n)                               \
    mpz_get_d(n)

// arithmetic
#define bn_add_to(n, ui)                               \
    mpz_add_ui(n, n, ui)

#define bn_incr(n)                               \
    mpz_add_ui(n, n, 1)

#define bn_add(n, m, o)                               \
    mpz_add(n, m, o)

#define bn_sub_from(n, ui)                               \
    mpz_sub_ui(n, n, ui)

#define bn_sub(n, m, o)                               \
    mpz_sub(n, m, o)

#define bn_mul_by(n, i)                               \
    mpz_mul_si(n, n, i)

#define bn_mul(n, m, o)                               \
    mpz_mul(n, m, o)

#define bn_div_by(n, ui)                               \
    mpz_fdiv_ui(n, n, ui)

#define bn_div(n, m, o)                               \
    mpz_fdiv_q(n, m, o)

// remainder from a long (returned as a long)
#define bn_reml(n, ui)                               \
    mpz_fdiv_ui(n, ui)

// remainder from a bignum
#define bn_rem(n, m, o)                               \
    mpz_mod(n, m, o)

// remainder from a bignum
#define bn_pow(n, m, ui)                               \
    mpz_pow_ui(n, m, ui)

// sign manipulation
#define bn_negate(n)                               \
    mpz_neg(n, n)

#define bn_set_neg(n, m)                               \
    mpz_neg(n, m)

#define bn_set_abs(n, m)                               \
    mpz_abs(n, m)

// returns the sign of n as -1, 0, or +1
#define bn_sgn(n)                               \
    mpz_sgn(n)

// comparisons
#define bn_cmp(n, m)                               \
    mpz_cmp(n, m)

#define bn_cmpl(n, i)                               \
    mpz_cmp_si(n, i)

#define bn_cmpd(n, d)                               \
    mpz_cmp_d(n, d)

// bit logic
#define bn_and(n, m)                               \
    mpz_and(n, m)

#define bn_or(n, m)                               \
    mpz_ior(n, m)

#define bn_xor(n, m)                               \
    mpz_xor(n, m)

// bit manipulation
#define bn_bitget(n, index)                               \
    mpz_tstbit(n, index)

#define bn_bitset(n, index, b)                               \
    {if (b) { mpz_setbit(n, index) } else { mpz_clrbit(n, index) }}

#define bn_bitflip(n, index)                               \
    mpz_combit(n, index)

// shifts
#define bn_sleft_by(n, ui)                               \
    mpz_mul_2exp(n, n, ui)

#define bn_set_sleft(n, m, ui)                               \
    mpz_mul_2exp(n, m, ui)

#define bn_sright_by(n, ui)                               \
    mpz_fdiv_q_2exp(n, n, ui)

#define bn_set_sright(n, m, ui)                               \
    mpz_fdiv_q_2exp(n, m, ui)

// printing
#define bn_print(n)                               \
    mpz_out_str(stdout, 10, n)
#define bn_print_hex(n)                               \
    mpz_out_str(stdout, 16, n)
#define bn_print_bin(n)                               \
    mpz_out_str(stdout, 2, n)

#define bn_fprint(stream, n)                               \
    mpz_out_str(stream, 10, n)
#define bn_fprint_hex(stream, n)                               \
    mpz_out_str(stream, 16, n)
#define bn_fprint_bin(stream, n)                               \
    mpz_out_str(stream, 2, n)

// size
#define bn_fits_in_long(n)                              \
    mpz_fits_slong_p(n) 

#endif