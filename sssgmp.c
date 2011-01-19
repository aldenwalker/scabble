#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "sssgmp.h"


void mpq_init(mpq_ptr a) {
  a->n->z = 0;
  a->d->z = 1;
}

void mpq_clear(mpq_ptr a) {
  //do nothing?
}

void mpq_out_str(FILE* fp, int base, mpq_ptr a) {
  if (fp == NULL) {
    if (a->d->z == 1L) {
      printf("%ld", a->n->z);
    } else {
      printf("%ld/%ld", a->n->z, a->d->z);
    }
  } else {
    if (a->d->z == 1L) {
      fprintf(fp, "%ld", a->n->z);
    } else {
      fprintf(fp, "%ld/%ld", a->n->z, a->d->z);
    }
  }
}

void mpq_set_si(mpq_ptr a, long b, long c) {
  a->n->z = b;
  a->d->z = c;
}

double mpq_get_d(mpq_ptr a) {
  return (double)a->n->z/(double)a->d->z;
}

void mpq_set(mpq_ptr a, mpq_ptr b) {
  a->n->z = b->n->z;
  a->d->z = b->d->z;
}

void mpq_neg(mpq_ptr a, mpq_ptr b) {
  a->d->z = b->d->z;
  a->n->z = -b->n->z;
}

int mpq_cmp_si(mpq_ptr a, long b, long c) {
  long right = a->d->z * b;
  long left = a->n->z * c;
  return (int)(left - right);
}

long mpq_cmp(mpq_ptr a, mpq_ptr b) {
  long right = a->d->z * b->n->z;
  long left = a->n->z * b->d->z;
  return (long)(left - right);
} 

mpz_ptr mpq_numref(mpq_ptr a) {
  return a->n;
}

mpz_ptr mpq_denref(mpq_ptr a) {
  return a->d;
}

void mpq_canonicalize(mpq_ptr a) {
  mpz_t gcd;
  //printf("canonicalizing %ld/%ld\n", a->n->z, a->d->z);
  mpz_gcd(gcd, mpq_numref(a), mpq_denref(a));
  a->n->z /= gcd->z;
  a->d->z /= gcd->z;
  if (a->d->z < 0) {
    a->n->z *= -1;
    a->d->z *= -1;
  }
  if (a->n->z == 0) {
    a->d->z = 1;
  }
  //printf("made it %ld/%ld\n", a->n->z, a->d->z);
}  
  
void mpq_mul(mpq_ptr a, mpq_ptr b, mpq_ptr c) {
  long tn1 = b->n->z;
  long tn2 = c->n->z;
  long td1 = b->d->z;
  long td2 = c->d->z;
  a->n->z = tn1 * tn2;
  a->d->z = td1 * td2;
  //printf("multiplying %ld/%ld by %ld/%ld and got %ld/%ld\n", tn1, td1, tn2, td2, a->n->z, a->d->z);
  mpq_canonicalize(a);
}

void mpq_add(mpq_ptr a, mpq_ptr b, mpq_ptr c) {
  long tn1 = b->n->z;
  long tn2 = c->n->z;
  long td1 = b->d->z;
  long td2 = c->d->z;  
  a->n->z = (td2*tn1) + (td1*tn2);
  a->d->z = td1 * td2;
  mpq_canonicalize(a);
}

void mpq_sub(mpq_ptr a, mpq_ptr b, mpq_ptr c) {
  long tn1 = b->n->z;
  long tn2 = c->n->z;
  long td1 = b->d->z;
  long td2 = c->d->z; 
  a->n->z = (td2*tn1) - (td1*tn2);
  a->d->z = td1 * td2; 
  mpq_canonicalize(a);
}

void mpq_div(mpq_ptr a, mpq_ptr b, mpq_ptr c) {
  long tn1 = b->n->z;
  long tn2 = c->n->z;
  long td1 = b->d->z;
  long td2 = c->d->z; 
  a->n->z = tn1*td2;
  a->d->z = td1*tn2;  
  //printf("dividing %ld/%ld by %ld/%ld and got %ld/%ld\n", b->n->z, b->d->z, c->n->z, c->d->z, a->n->z, a->d->z);
  mpq_canonicalize(a);
}

void mpq_swap(mpq_ptr a, mpq_ptr b) {
  long tn = a->n->z;
  long td = a->d->z;
  a->n->z = b->n->z;
  a->d->z = b->d->z;
  b->n->z = tn;
  b->d->z = td;
}

void mpq_set_z(mpq_ptr a, mpz_ptr b) {
  a->n->z = b->z;
  a->d->z = 1L;
}







void mpz_init(mpz_ptr a) {
  //do nothing?
}

void mpz_clear(mpz_ptr a) {
  //do nothing?
}

void mpz_set(mpz_ptr a, mpz_ptr b) {
  a->z = b->z;
}

void mpz_set_si(mpz_ptr a, long b) {
  a->z = b;
}

void mpz_lcm(mpz_ptr a, mpz_ptr b, mpz_ptr c) {
  mpz_t gcd;
  mpz_gcd(gcd, b, c);
  a->z = (b->z * c->z)/(gcd->z);
}

void mpz_gcd(mpz_ptr a, mpz_ptr b, mpz_ptr c) {
  long t1 = (b->z > c->z ? b->z : c->z); //max
  long t2 = (b->z > c->z ? c->z : b->z); //min
  long t3;
  while (t2 != 0) {
    t3 = t1 % t2;
    t1 = t2;
    t2 = t3;
  }
  a->z = t1;
}
  
long mpz_get_si(mpz_ptr a) {
  return a->z;
}

long mpz_cmp_si(mpz_ptr a, long b) {
  return a->z - b;
}
void mpz_fdiv_q(mpz_ptr a, mpz_ptr b, mpz_ptr c) {
  if (b->z == 0) {
    a->z = 0;
    return;
  }
  if ((b->z > 0 && c->z > 0) || (b->z < 0 && c->z < 0)) {  //regular is correct
    a->z = (b->z)/(c->z);
  } else {
    if ((b->z)%(c->z) == 0) {
      a->z = (b->z)/(c->z);
    } else {
      a->z = (b->z)/(c->z) - 1;
    }
  }
}


int cmp_doubles(double a, double b) {
  if (fabs(a-b) < 1e-10) {
    return 1;
  } else {
    return 0;
  }
}

double round_double(double a) {
  if (a >= 0) {
    if (fmod(a, 1) > 0.5) {
      return ceil(a);
    } else {
      return floor(a);
    }
  } else {
    if (fabs(fmod(a,1)) > 0.5) {
      return floor(a);
    } else {
      return ceil(a);
    }
  } 
}

double sanize_double(double a) {
  if (cmp_doubles(round_double(a), a)) {
    return round_double(a);
  } else {
    return a;
  }
}








