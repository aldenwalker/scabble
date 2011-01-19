#ifndef _sss_gmp_defined
#define _sss_gmp_defined

/*
these functions are meant to be limited precision copies of the GMP functions
so that we could make it unlimited precision just by changing the #include's
*/

#include <stdio.h>

typedef struct {
  long z;
} __mpz_struct;

typedef __mpz_struct mpz_t[1];
typedef __mpz_struct* mpz_ptr;

typedef struct {
  mpz_t n;
  mpz_t d;
} __mpq_struct;

typedef __mpq_struct mpq_t[1];
typedef __mpq_struct* mpq_ptr;

void mpq_init(mpq_ptr a);
void mpq_clear(mpq_ptr a);
void mpq_out_str(FILE* fp, int base, mpq_ptr a);
void mpq_set_si(mpq_ptr a, long b, long c);

double mpq_get_d(mpq_ptr a);
void mpq_set(mpq_ptr a, mpq_ptr b);
void mpq_neg(mpq_ptr a, mpq_ptr b);
int mpq_cmp_si(mpq_ptr a, long b, long c);
long mpq_cmp(mpq_ptr a, mpq_ptr b);
mpz_ptr mpq_numref(mpq_ptr a);
mpz_ptr mpq_denref(mpq_ptr a);
void mpq_canonicalize(mpq_ptr a);
void mpq_mul(mpq_ptr a, mpq_ptr b, mpq_ptr c);
void mpq_add(mpq_ptr a, mpq_ptr b, mpq_ptr c);
void mpq_sub(mpq_ptr a, mpq_ptr b, mpq_ptr c);
void mpq_div(mpq_ptr a, mpq_ptr b, mpq_ptr c);
void mpq_swap(mpq_ptr a, mpq_ptr b);
void mpq_set_z(mpq_ptr a, mpz_ptr b);

void mpz_init(mpz_ptr a);
void mpz_clear(mpz_ptr a);
void mpz_set(mpz_ptr a, mpz_ptr b);
void mpz_set_si(mpz_ptr a, long b);
void mpz_lcm(mpz_ptr a, mpz_ptr b, mpz_ptr c);
void mpz_gcd(mpz_ptr a, mpz_ptr b, mpz_ptr c);
long mpz_get_si(mpz_ptr a);
long mpz_cmp_si(mpz_ptr a, long b);
void mpz_fdiv_q(mpz_ptr a, mpz_ptr b, mpz_ptr c);

int cmp_doubles(double a, double b);
double round_double(double a);
double sanize_double(double a);




#endif
