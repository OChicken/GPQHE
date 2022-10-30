#include "ckks.h"
//#include "controller.h"
//#include "params.h"
//#include "precomp.h"

#include <fcntl.h>   /* open & creat, return int */
#include <unistd.h>  /* read & write, return ssize_t */

/* BLAS */
#include <cblas.h>
#include <lapacke.h>

#include <pmu.h>

#ifndef GPQHE_PRECOMP
#define GPQHE_PRECOMP 0
#endif

#define NTESTS 2

/*
 * #pragma cling load("libgcrypt.so");
 * #pragma cling load("libckks.so");
 * #include <gcrypt.h>
 * #include "src/ckks.h"
*/

__attribute__((unused)) static void
test_convert(void)
{
  size_t l=15;
  unsigned char *buf = (unsigned char *)malloc((l/8+1)*sizeof(unsigned char));
  buf[0]='A';
  buf[1]='B';
  printf("l=%ld: %li\n", l  , ntdk_buf_to_u64(buf, l));
  printf("l=%ld: %li\n", l-1, ntdk_buf_to_u64(buf, l-1));
  l = 5*8;
  buf = (unsigned char *)malloc((l/8+1)*sizeof(unsigned char));
  buf[0]='c',buf[1]='h',buf[2]='a',buf[3]='r',buf[4]='\a';
  printf("buf=\"char\\a\"=%ld\n",ntdk_buf_to_u64(buf,l));
  size_t num = (1<< 6) + (1<< 5) + (1<< 1) + (1<<0) + \
               (1<<14) + (1<<13) + (1<<11) + \
               (1<<22) + (1<<21) + (1<<16) + \
               (1<<30) + (1<<29) + (1<<28) + (1<<25) + \
               ((uint64_t)1<<34) + ((uint64_t)1<<33) + ((uint64_t)1<<32);
  assert(ntdk_buf_to_u64(buf,l)==num);
  free(buf);
  printf("\033[1m%s pass.\033[0m\n", __func__);
}

__attribute__((unused)) static void
test_bit_operations(void)
{
  size_t size;
  size = 4;
  long *c = (long*) calloc(size, sizeof(long));
  long *r = (long*) calloc(size, sizeof(long));
  long *a = (long*) calloc(size, sizeof(long));
  c[0]=0,c[1]=1,c[2]=2,c[3]=3;
  a[0]=0,a[1]=2,a[2]=1,a[3]=3;
  for (size_t i=0; i<size; i++)
    r[i] = c[mpi_br_width(i, ntdk_log2(size))];
  assert(blas_lvec_equal(size, a, r)==1);
  free(c);
  free(r);
  free(a);
  size = 8;
  c = (long*) calloc(size, sizeof(long)); /* coeffs */
  r = (long*) calloc(size, sizeof(long)); /* return */
  a = (long*) calloc(size, sizeof(long)); /* answer */
  c[0]=0,c[1]=1,c[2]=2,c[3]=3;c[4]=4,c[5]=5,c[6]=6,c[7]=7;
  a[0]=0,a[1]=4,a[2]=2,a[3]=6;a[4]=1,a[5]=5,a[6]=3,a[7]=7;
  for (size_t i=0; i<size; i++)
    r[i] = c[mpi_br_width(i, ntdk_log2(size))];
  assert(blas_lvec_equal(size, a, r)==1);
  c[0]=1,c[1]=3,c[2]=5,c[3]= 7;c[4]=9,c[5]=11,c[6]=13,c[7]=15;
  a[0]=1,a[1]=9,a[2]=5,a[3]=13;a[4]=3,a[5]=11,a[6]= 7,a[7]=15;
  for (size_t i=0; i<size; i++)
    r[i] = c[mpi_br_width(i, ntdk_log2(size))];
  assert(blas_lvec_equal(size, a, r)==1);
  free(c);
  free(r);
  free(a);
  size = 16;
  c = (long*) calloc(size, sizeof(long));
  r = (long*) calloc(size, sizeof(long));
  a = (long*) calloc(size, sizeof(long));
  for (size_t i=0; i<size; i++)
    c[i]=i;
  a[0]=0,a[1]=8,a[ 2]=4,a[ 3]=12,a[ 4]=2,a[ 5]=10,a[ 6]=6,a[ 7]=14;
  a[8]=1,a[9]=9,a[10]=5,a[11]=13,a[12]=3,a[13]=11,a[14]=7,a[15]=15;
  for (size_t i=0; i<size; i++)
    r[i] = c[mpi_br_width(i, ntdk_log2(size))];
  assert(blas_lvec_equal(size, a, r)==1);
  free(c);
  free(r);
  free(a);
  printf("\033[1m%s pass.\033[0m\n", __func__);
}

__attribute__((unused)) static void
test_crt(void)
{
  crt_ctx crt;
  size_t num_primes = 4;
  size_t size_p = 48;
  crt_build(&crt, size_p, num_primes);
  gcry_mpi_t original = gcry_mpi_new(0);
  gcry_mpi_t reverse = gcry_mpi_new(0);
  gcry_mpi_t *vals = gcry_calloc(num_primes, sizeof(gcry_mpi_t));
  for (size_t i=0; i<num_primes; i++)
    vals[i] = gcry_mpi_new(0);
  original = gcry_int64_to_mpi(178);
  crt_fwd(vals, original, &crt);
  crt_inv(&reverse, vals, &crt);
  for (size_t i=0; i<num_primes; i++)
    printf("%lu ", gcry_mpi_to_int64(vals[i]));
  printf("\n");
  printf("reverse=%ld\n", gcry_mpi_to_int64(reverse));
  assert(gcry_mpi_to_int64(original) == gcry_mpi_to_int64(reverse));
  free(vals);
  printf("\033[1m%s pass.\033[0m\n", __func__);
}

__attribute__((unused)) static void
test_ntt(void)
{
  ntt_ctx ntt;
  gcry_mpi_t modulus = gcry_mpi_set_ui(NULL, 73);
  gcry_mpi_t g = gcry_mpi_set_ui(NULL, ntdk_m_root_of_unity(72, 73));
  ntt_build(&ntt, 4, modulus, g);
  printf("ntt.zetas: ");
  for (size_t i=0; i<4; i++)
    printf("%ld ", gcry_mpi_to_int64(ntt.zetas[i]));
  printf("\n");
  size_t m = 4;
  gcry_mpi_t *x    = gcry_calloc(m, sizeof(gcry_mpi_t));
  gcry_mpi_t *X    = gcry_calloc(m, sizeof(gcry_mpi_t));
  gcry_mpi_t *inv  = gcry_calloc(m, sizeof(gcry_mpi_t));
  for (size_t i=0; i<m; i++)
  {
    x[i] = gcry_mpi_new(0);
    X[i] = gcry_mpi_new(0);
    inv[i] = gcry_mpi_new(0);
  }
  /* unsigned case */
  x[0] = gcry_int64_to_mpi(0);
  x[1] = gcry_int64_to_mpi(1);
  x[2] = gcry_int64_to_mpi(4);
  x[3] = gcry_int64_to_mpi(5);
  ntt_fwd(m, X, x, modulus, &ntt);
  ntt_inv(m, inv, X, modulus, &ntt);
  printf("[0,1,4,+5]'s X: ");
  for (size_t i=0; i<m; i++)
    printf("%ld ", gcry_mpi_to_int64(X[i]));
  printf("\n");
  assert(blas_mpivec_equal(m, inv, x)==1);
  /* signed case */
  x[0] = gcry_int64_to_mpi(0);
  x[1] = gcry_int64_to_mpi(1);
  x[2] = gcry_int64_to_mpi(4);
  x[3] = gcry_int64_to_mpi(-5);
  ntt_fwd(m, X, x, modulus, &ntt);
  ntt_inv(m, inv, X, modulus, &ntt);
  for (size_t i=0; i<m; i++)
    ntdk_smod(&inv[i], modulus);
  printf("[0,1,4,-5]'s X: ");
  for (size_t i=0; i<m; i++)
    printf("%ld ", gcry_mpi_to_int64(X[i]));
  printf("\n");
  assert(blas_mpivec_equal(m, inv, x)==1);
  gcry_free(x);
  gcry_free(X);
  gcry_free(inv);
  printf("\033[1m%s pass.\033[0m\n", __func__);
}

__attribute__((unused)) static void
test_fft(void)
{
  size_t degree = 16;
  size_t num_slots = degree/2;
  fft_ctx fft;
  fft_build(&fft, degree);
  int m=4;
  COMPLEX *coeffs = (COMPLEX *)calloc(m, sizeof(COMPLEX));
  COMPLEX *fwd    = (COMPLEX *)calloc(m, sizeof(COMPLEX));
  COMPLEX *inv    = (COMPLEX *)calloc(m, sizeof(COMPLEX));
  COMPLEX *answer = (COMPLEX *)calloc(m, sizeof(COMPLEX));
  coeffs[0]=0,coeffs[1]=1,coeffs[2]=4,coeffs[3]=5;
  answer[0]=10,answer[1]=-4-4*I,answer[2]=-2,answer[3]=-4+4*I;
  fft_fwd(m, fwd, coeffs, &fft);
  fft_inv(m, inv, fwd, &fft);
  printf("||fwd-answer|| = %g\n", blas_dznrmmax_dist(m, fwd, answer));
  printf("||inv-coeffs|| = %g\n", blas_dznrmmax_dist(m, inv, coeffs));
  assert(blas_zvec_equal(m, fwd, answer, 1e-15)==1);
  assert(blas_zvec_equal(m, inv, coeffs, __DBL_EPSILON__)==1);
  free(coeffs);
  free(fwd);
  free(inv);
  free(answer);
  coeffs = (COMPLEX *)calloc(num_slots, sizeof(COMPLEX));
  fwd    = (COMPLEX *)calloc(num_slots, sizeof(COMPLEX));
  inv    = (COMPLEX *)calloc(num_slots, sizeof(COMPLEX));
  sample_z01vec(num_slots, coeffs);
  fft_fwd(num_slots, fwd, coeffs, &fft);
  fft_inv(num_slots, inv, fwd, &fft);
  printf("||inv-coeffs|| = %g\n", blas_dznrmmax_dist(m, inv, coeffs));
  printf("\033[1m%s pass.\033[0m\n", __func__);
}

__attribute__((unused)) static void
test_embedding(void)
{
  size_t degree = 16;
  size_t m = degree/2;
  fft_ctx fft;
  fft_build(&fft, degree);
  COMPLEX *coeffs = (COMPLEX *)calloc(m, sizeof(COMPLEX));
  COMPLEX *emb    = (COMPLEX *)calloc(m, sizeof(COMPLEX));
  COMPLEX *inv    = (COMPLEX *)calloc(m, sizeof(COMPLEX));
  COMPLEX *eval   = (COMPLEX *)calloc(m, sizeof(COMPLEX));
  coeffs[0]=10,coeffs[1]=34,coeffs[2]=71,coeffs[3]=31;
  coeffs[4]=1, coeffs[5]=2, coeffs[6]=3, coeffs[7]=4;
  emb_fwd(m, emb, coeffs, &fft);
  emb_inv(m, inv, emb, &fft);
  printf("||inv-coeffs|| = %g\n", blas_dznrmmax_dist(m, inv, coeffs));
  assert(blas_zvec_equal(m, inv, coeffs, 1e-13)==1);
  printf("emb_fwd and emb_inv are dual operations.\n");
  size_t power=1,idx=0;
  for (size_t i=1; i<fft.L; i+=4)
  {
    double angle = 2*M_PI*power/fft.L;
    COMPLEX zeta = cos(angle) + I*sin(angle);
    COMPLEX ret = coeffs[m-1];
    for (int j=m-2; j>-1; j--)
      ret = ret*zeta + coeffs[j];
    eval[idx++] = ret;
    power = (power*5)%fft.L;
  }
  printf("||emb-eval|| = %g\n", blas_dznrmmax_dist(m, emb, eval));
  assert(blas_zvec_equal(m, emb, eval, 1e-13));
  printf("embedding matches the evaluations of the roots of unity at indices that are 1mod4.\n");
  printf("\033[1m%s pass.\033[0m\n", __func__);
}

#if 0
void test_poly_mod()
{
  size_t degree = 5;
  gcry_mpi_t modulus = gcry_mpi_set_ui(NULL, 60);
  gcry_mpi_t *poly0  = (gcry_mpi_t *)calloc(degree, sizeof(gcry_mpi_t));
  gcry_mpi_t *poly1  = (gcry_mpi_t *)calloc(degree, sizeof(gcry_mpi_t));
  gcry_mpi_t *poly2  = (gcry_mpi_t *)calloc(degree, sizeof(gcry_mpi_t));
  gcry_mpi_t *answer = (gcry_mpi_t *)calloc(degree, sizeof(gcry_mpi_t));
  for (size_t i=0; i<degree; i++)
  {
    poly0[i] = gcry_mpi_new(0);
    poly1[i] = gcry_mpi_new(0);
    poly2[i] = gcry_mpi_new(0);
    answer[i] = gcry_mpi_new(0);
  }
  poly1[0]=gcry_int64_to_mpi(57)   , answer[0]=gcry_int64_to_mpi(57);
  poly1[1]=gcry_int64_to_mpi(-34)  , answer[1]=gcry_int64_to_mpi(26);
  poly1[2]=gcry_int64_to_mpi(100)  , answer[2]=gcry_int64_to_mpi(40);
  poly1[3]=gcry_int64_to_mpi(1000) , answer[3]=gcry_int64_to_mpi(40);
  poly1[4]=gcry_int64_to_mpi(-7999), answer[4]=gcry_int64_to_mpi(41);
  poly_mod(degree, poly1, modulus);
  assert(blas_mpivec_equal(degree, poly1, answer));
  printf("\033[1m%s pass.\033[0m\n", __func__);
}

void test_poly_add()
{
  size_t degree = 5;
  gcry_mpi_t modulus = gcry_mpi_set_ui(NULL, 60);
  gcry_mpi_t *poly0  = (gcry_mpi_t *)calloc(degree, sizeof(gcry_mpi_t));
  gcry_mpi_t *poly1  = (gcry_mpi_t *)calloc(degree, sizeof(gcry_mpi_t));
  gcry_mpi_t *poly2  = (gcry_mpi_t *)calloc(degree, sizeof(gcry_mpi_t));
  gcry_mpi_t *answer = (gcry_mpi_t *)calloc(degree, sizeof(gcry_mpi_t));
  for (size_t i=0; i<degree; i++)
  {
    poly0[i] = gcry_mpi_new(0);
    poly1[i] = gcry_mpi_new(0);
    poly2[i] = gcry_mpi_new(0);
    answer[i] = gcry_mpi_new(0);
  }
  /* test add and mod */
  poly1[0]=gcry_int64_to_mpi(0) , poly2[0]=gcry_int64_to_mpi(1), answer[0]=gcry_int64_to_mpi(1);
  poly1[1]=gcry_int64_to_mpi(1) , poly2[1]=gcry_int64_to_mpi(2), answer[1]=gcry_int64_to_mpi(3);
  poly1[2]=gcry_int64_to_mpi(4) , poly2[2]=gcry_int64_to_mpi(4), answer[2]=gcry_int64_to_mpi(8);
  poly1[3]=gcry_int64_to_mpi(5) , poly2[3]=gcry_int64_to_mpi(3), answer[3]=gcry_int64_to_mpi(8);
  poly1[4]=gcry_int64_to_mpi(59), poly2[4]=gcry_int64_to_mpi(2), answer[4]=gcry_int64_to_mpi(1);
  poly_xpy(degree, poly0, poly1, poly2, modulus, 0);
  assert(blas_mpivec_equal(degree, poly0, answer));
  printf("\033[1m%s pass.\033[0m\n", __func__);
}

void test_poly_sub()
{
  size_t degree = 5;
  gcry_mpi_t modulus = gcry_mpi_set_ui(NULL, 60);
  gcry_mpi_t *poly0  = (gcry_mpi_t *)calloc(degree, sizeof(gcry_mpi_t));
  gcry_mpi_t *poly1  = (gcry_mpi_t *)calloc(degree, sizeof(gcry_mpi_t));
  gcry_mpi_t *poly2  = (gcry_mpi_t *)calloc(degree, sizeof(gcry_mpi_t));
  gcry_mpi_t *answer = (gcry_mpi_t *)calloc(degree, sizeof(gcry_mpi_t));
  for (size_t i=0; i<degree; i++)
  {
    poly0[i] = gcry_mpi_new(0);
    poly1[i] = gcry_mpi_new(0);
    poly2[i] = gcry_mpi_new(0);
    answer[i] = gcry_mpi_new(0);
  }
  /* test sub and mod */
  poly1[0]=gcry_int64_to_mpi(0) , poly2[0]=gcry_int64_to_mpi(1), answer[0]=gcry_int64_to_mpi(59);
  poly1[1]=gcry_int64_to_mpi(1) , poly2[1]=gcry_int64_to_mpi(2), answer[1]=gcry_int64_to_mpi(59);
  poly1[2]=gcry_int64_to_mpi(4) , poly2[2]=gcry_int64_to_mpi(4), answer[2]=gcry_int64_to_mpi(0);
  poly1[3]=gcry_int64_to_mpi(5) , poly2[3]=gcry_int64_to_mpi(3), answer[3]=gcry_int64_to_mpi(2);
  poly1[4]=gcry_int64_to_mpi(59), poly2[4]=gcry_int64_to_mpi(2), answer[4]=gcry_int64_to_mpi(57);
  poly_xsy(degree, poly0, poly1, poly2, modulus, 0);
  assert(blas_mpivec_equal(degree, poly0, answer));
  printf("\033[1m%s pass.\033[0m\n", __func__);
}

void test_poly_scal()
{
  size_t degree = 5;
  gcry_mpi_t modulus = gcry_mpi_set_ui(NULL, 60);
  gcry_mpi_t *poly1  = gcry_calloc(degree, sizeof(gcry_mpi_t));
  gcry_mpi_t *answer = gcry_calloc(degree, sizeof(gcry_mpi_t));
  for (size_t i=0; i<degree; i++)
  {
    poly1[i] = gcry_mpi_new(0);
    answer[i] = gcry_mpi_new(0);
  }
  /* test scal and mod */
  poly1[0]=gcry_int64_to_mpi(0) , answer[0]=gcry_int64_to_mpi(0);
  poly1[1]=gcry_int64_to_mpi(1) , answer[1]=gcry_int64_to_mpi(59);
  poly1[2]=gcry_int64_to_mpi(4) , answer[2]=gcry_int64_to_mpi(56);
  poly1[3]=gcry_int64_to_mpi(5) , answer[3]=gcry_int64_to_mpi(55);
  poly1[4]=gcry_int64_to_mpi(59), answer[4]=gcry_int64_to_mpi(1);
  poly_scal(degree, poly1, -1);
  poly_mod(degree, poly1, modulus);
  for (size_t i=0; i<degree; i++)
    printf("%ld ", gcry_mpi_to_int64(poly1[i]));
  printf("\n");
  for (size_t i=0; i<degree; i++)
    printf("%ld ", gcry_mpi_to_int64(answer[i]));
  printf("\n");
  if (gcry_mpi_cmp_ui(poly1[0], 0)==0)
    printf("poly1 correct\n");
  if (gcry_mpi_cmp_ui(answer[0], 0)==0)
    printf("answer correct\n");
  if (gcry_mpi_cmp(poly1[0], answer[0])==0)
    printf("poly1[0]==answer[0]\n");
  printf("\033[1m%s pass.\033[0m\n", __func__);
}

void test_poly_rot()
{
  /* poly_rot */
  size_t degree = 4;
  gcry_mpi_t *poly1 = gcry_calloc(degree, sizeof(gcry_mpi_t));
  gcry_mpi_t *poly2 = gcry_calloc(degree, sizeof(gcry_mpi_t));
  gcry_mpi_t *answer = gcry_calloc(degree, sizeof(gcry_mpi_t));
  for (size_t i=0; i<degree; i++)
  {
    poly1[i] = gcry_mpi_new(0);
    poly2[i] = gcry_mpi_new(0);
    answer[i] = gcry_mpi_new(0);
  }
  poly1[0]=gcry_int64_to_mpi(0), answer[0]=gcry_int64_to_mpi(0);
  poly1[1]=gcry_int64_to_mpi(1), answer[1]=gcry_int64_to_mpi(-1);
  poly1[2]=gcry_int64_to_mpi(4), answer[2]=gcry_int64_to_mpi(4);
  poly1[3]=gcry_int64_to_mpi(59),answer[3]=gcry_int64_to_mpi(-59);
  poly_rot(degree, poly2, poly1, 3);
  assert(blas_mpivec_equal(degree, poly2, answer));
  printf("\033[1m%s pass.\033[0m\n", __func__);
}

void test_polynomial()
{
  size_t degree = 5;
  gcry_mpi_t modulus = gcry_mpi_set_ui(NULL, 60);
  gcry_mpi_t *poly0  = (gcry_mpi_t *)calloc(degree, sizeof(gcry_mpi_t));
  gcry_mpi_t *poly1  = (gcry_mpi_t *)calloc(degree, sizeof(gcry_mpi_t));
  gcry_mpi_t *poly2  = (gcry_mpi_t *)calloc(degree, sizeof(gcry_mpi_t));
  gcry_mpi_t *answer = (gcry_mpi_t *)calloc(degree, sizeof(gcry_mpi_t));
  for (size_t i=0; i<degree; i++)
  {
    poly0[i] = gcry_mpi_new(0);
    poly1[i] = gcry_mpi_new(0);
    poly2[i] = gcry_mpi_new(0);
    answer[i] = gcry_mpi_new(0);
  }
  polyd polyd1;
  #if 0
  /* round and mod */
  degree = 5;
  poly1.d = degree;
  polyd1.d = degree;
  poly1.coeffs = (int64_t *)calloc(degree, sizeof(int64_t));
  poly2.coeffs = (int64_t *)calloc(degree, sizeof(int64_t));
  polyd1.coeffs = (double *)calloc(degree, sizeof(double));
  coeffs = (int64_t *)calloc(degree, sizeof(int64_t));
  polyd1.coeffs[0]=0.51,polyd1.coeffs[1]=-3.2,polyd1.coeffs[2]=54.666,polyd1.coeffs[3]=39.01,polyd1.coeffs[4]=0;
  coeffs[0]=1,coeffs[1]=-3,coeffs[2]=55,coeffs[3]=39,coeffs[4]=0;
  poly_round(&poly1, &polyd1);
  assert(blas_lvec_equal(degree, poly1.coeffs, coeffs));
  poly1.coeffs[0]=57,poly1.coeffs[1]=-34,poly1.coeffs[2]=100,poly1.coeffs[3]=1000,poly1.coeffs[4]=-7999;
  coeffs[0]=57,coeffs[1]=26,coeffs[2]=40,coeffs[3]=40,coeffs[4]=41;
  poly_mod(&poly1, modulus);
  assert(blas_lvec_equal(degree, poly1.coeffs, coeffs));
  printf("test poly round and mod passed.\n");
  #endif

  /* poly_eval */
  degree = 5;
  poly1 = gcry_calloc(degree, sizeof(gcry_mpi_t));
  for (size_t i=0; i<degree; i++)
  {
    poly1[i] = gcry_mpi_new(0);
    poly1[i] = gcry_int64_to_mpi(i);
  }
  assert(gcry_mpi_to_int64(poly_eval(degree, poly1, 3))==426);
  printf(("test poly evaluation passed.\n"));
  printf("\033[1m%s pass.\033[0m\n", __func__);
}

void test_poly_multiply()
{
  size_t degree = 4;
  gcry_mpi_t modulus = gcry_mpi_set_ui(NULL, 60);
  gcry_mpi_t *poly0  = gcry_calloc(degree, sizeof(gcry_mpi_t));
  gcry_mpi_t *poly1  = gcry_calloc(degree, sizeof(gcry_mpi_t));
  gcry_mpi_t *poly2  = gcry_calloc(degree, sizeof(gcry_mpi_t));
  gcry_mpi_t *answer = gcry_calloc(degree, sizeof(gcry_mpi_t));
  for (size_t i=0; i<degree; i++)
  {
    poly0[i] = gcry_mpi_new(0);
    poly1[i] = gcry_mpi_new(0);
    poly2[i] = gcry_mpi_new(0);
    answer[i] = gcry_mpi_new(0);
  }
  /* test multiply naive, and then mod */
  printf("test poly_xmy\n");
  poly1[0]=gcry_int64_to_mpi(0) , poly2[0]=gcry_int64_to_mpi(1);
  poly1[1]=gcry_int64_to_mpi(1) , poly2[1]=gcry_int64_to_mpi(2);
  poly1[2]=gcry_int64_to_mpi(4) , poly2[2]=gcry_int64_to_mpi(4);
  poly1[3]=gcry_int64_to_mpi(5) , poly2[3]=gcry_int64_to_mpi(3);

  /* poly_mod: take big mod 73 */
  modulus = gcry_mpi_set_ui(NULL, 73);
  answer[0]=gcry_int64_to_mpi(44);
  answer[1]=gcry_int64_to_mpi(42);
  answer[2]=gcry_int64_to_mpi(64);
  answer[3]=gcry_int64_to_mpi(17);
  poly_xmy(degree, poly0, poly1, poly2, modulus, 0);
  assert(blas_mpivec_equal(degree, poly0, answer));
  printf("  take big mod 73 passed.\n");

  /* poly_smod: take sighed mod */
  poly0  = gcry_calloc(degree, sizeof(gcry_mpi_t));
  for (size_t i=0; i<degree; i++)
    poly0[i] = gcry_mpi_new(0);
  poly_xmy(degree, poly0, poly1, poly2, modulus, 1);
  answer[0]=gcry_int64_to_mpi(-29);
  answer[1]=gcry_int64_to_mpi(-31);
  answer[2]=gcry_int64_to_mpi(-9);
  answer[3]=gcry_int64_to_mpi(17);
  assert(blas_mpivec_equal(degree, poly0, answer));
  printf("  take signed mod 73 passed.\n");

  /* poly_mod: take big mod 60 */
  degree = 5;
  poly0  = gcry_calloc(degree, sizeof(gcry_mpi_t));
  poly1  = gcry_calloc(degree, sizeof(gcry_mpi_t));
  poly2  = gcry_calloc(degree, sizeof(gcry_mpi_t));
  answer = gcry_calloc(degree, sizeof(gcry_mpi_t));
  for (size_t i=0; i<degree; i++)
  {
    poly0[i] = gcry_mpi_new(0);
    poly1[i] = gcry_mpi_new(0);
    poly2[i] = gcry_mpi_new(0);
    answer[i] = gcry_mpi_new(0);
  }
  modulus = gcry_mpi_set_ui(NULL, 60);
  poly1[0]=gcry_int64_to_mpi(0) , poly2[0]=gcry_int64_to_mpi(1), answer[0]=gcry_int64_to_mpi(28);
  poly1[1]=gcry_int64_to_mpi(1) , poly2[1]=gcry_int64_to_mpi(2), answer[1]=gcry_int64_to_mpi(42);
  poly1[2]=gcry_int64_to_mpi(4) , poly2[2]=gcry_int64_to_mpi(4), answer[2]=gcry_int64_to_mpi(59);
  poly1[3]=gcry_int64_to_mpi(5) , poly2[3]=gcry_int64_to_mpi(3), answer[3]=gcry_int64_to_mpi(19);
  poly1[4]=gcry_int64_to_mpi(59), poly2[4]=gcry_int64_to_mpi(2), answer[4]=gcry_int64_to_mpi(28);
  poly_xmy(degree, poly0, poly1, poly2, modulus, 0);
  assert(blas_mpivec_equal(degree, poly0, answer));
  printf("  take big mod 60 passed.\n");
  printf("test poly_xmy x*y (naive) and mod & smod passed.\n");

  /* test multiply crt */
  size_t pbnd = 59;
  size_t log_d = 2;
  size_t log_modulus = 100;
  degree = 1 << log_d;
  gcry_mpi_lshift(modulus, gcry_mpi_set_ui(NULL, 1), log_modulus);
  size_t num_primes = (2 + log_d + 4 * log_modulus + pbnd - 1) / pbnd;
  poly0  = gcry_calloc(degree, sizeof(gcry_mpi_t));
  poly1  = gcry_calloc(degree, sizeof(gcry_mpi_t));
  poly2  = gcry_calloc(degree, sizeof(gcry_mpi_t));
  answer = gcry_calloc(degree, sizeof(gcry_mpi_t));
  for (size_t i=0; i<degree; i++)
  {
    poly0[i] = gcry_mpi_new(0);
    poly1[i] = gcry_mpi_new(0);
    poly2[i] = gcry_mpi_new(0);
    answer[i] = gcry_mpi_new(0);
  }
  poly1[0]=gcry_int64_to_mpi(0) , poly2[0]=gcry_int64_to_mpi(1), answer[0]=gcry_int64_to_mpi(-29);
  poly1[1]=gcry_int64_to_mpi(1) , poly2[1]=gcry_int64_to_mpi(2), answer[1]=gcry_int64_to_mpi(-31);
  poly1[2]=gcry_int64_to_mpi(4) , poly2[2]=gcry_int64_to_mpi(4), answer[2]=gcry_int64_to_mpi(-9);
  poly1[3]=gcry_int64_to_mpi(5) , poly2[3]=gcry_int64_to_mpi(3), answer[3]=gcry_int64_to_mpi(17);
  /* the answer: use naive method */
  poly_xmy(degree, poly0, poly1, poly2, modulus, 1);
  assert(blas_mpivec_equal(degree, poly0, answer));
  /* use CRT */
  crt_ctx crt;
  crt_build(&crt, pbnd, num_primes);
  ntt_ctx *ntts = (ntt_ctx *)calloc(num_primes, sizeof(ntt_ctx));
  for (size_t i=0; i<num_primes; i++)
    ntt_build(&ntts[i], degree, crt.p[i], crt.g[i]);
  poly_xmy_crt(degree, poly0, poly1, poly2, &crt, ntts);
  poly_smod(degree, poly0, modulus);
  for (size_t i=0; i<degree; i++)
    printf("%ld ", gcry_mpi_to_int64(poly0[i]));
  printf("\n");
  printf("\033[1m%s pass.\033[0m\n", __func__);
}

#endif
__attribute__((unused)) static void
test_encode()
{
  TEST_BEGIN();
  size_t d = 16;
  size_t m = 4;
  double pt[GPQHE_INDCPA_MSGBYTES];
  COMPLEX msg[m];
  COMPLEX dcd[m];
  memcpy(msg, (COMPLEX []){1+I*2, 2+I*3, 3+I*4, 4+I*5}, m*sizeof(COMPLEX));
  printf("msg = \033[1m\033[4m[\033[0m ");
  for (size_t i=0; i<m; i++)
    printf("%g+%gj ", creal(msg[i]), cimag(msg[i]));
  printf("\033[1m\033[4m]\033[0m\n");
  he_ckks_encode(pt, m, msg);
  printf("pt  = \033[1m\033[4m[\033[0m ");
  for (size_t i=0; i<d; i++)
    printf("%.0f ", pt[i]);
  printf("\033[1m\033[4m]\033[0m\n");
  he_ckks_decode(m, dcd, pt);
  printf("dcd = \033[1m\033[4m[\033[0m ");
  for (size_t i=0; i<m; i++)
    printf("%g+%gj ", creal(dcd[i]), cimag(dcd[i]));
  printf("\033[1m\033[4m]\033[0m\n");
  printf("||msg - msg_dec|| = %g\n", blas_dznrmmax_dist(m, msg, dcd));
#if 0
  /* prod */
  msg = (COMPLEX *)calloc(m, sizeof(COMPLEX));
  dcd = (COMPLEX *)calloc(m, sizeof(COMPLEX));
  ckks_pt pt1,pt2;
  gcry_mpi_t *mpi_pt_prod = gcry_calloc(d, sizeof(gcry_mpi_t));
  gcry_mpi_t *mpi_pt1     = gcry_calloc(d, sizeof(gcry_mpi_t));
  gcry_mpi_t *mpi_pt2     = gcry_calloc(d, sizeof(gcry_mpi_t));
  COMPLEX *msg1   = (COMPLEX *)calloc(m, sizeof(COMPLEX));
  COMPLEX *msg2   = (COMPLEX *)calloc(m, sizeof(COMPLEX));
  msg1[0]=1+I*2, msg1[1]=2+I*3, msg1[2]=3+I*4, msg1[3]=4+I*5;
  msg2[0]=5+I*6, msg2[1]=6+I*7, msg2[2]=7+I*8, msg2[3]=8+I*9;
  blas_zxmy(m, msg, msg1, msg2);
  he_ckks_encode(&pt, m, msg, &ctx);

  he_ckks_encode(&pt1, m, msg1, &ctx);
  he_ckks_encode(&pt2, m, msg2, &ctx);
  for (size_t i=0; i<d; i++)
  {
    mpi_pt_prod[i] = gcry_mpi_new(0);
    mpi_pt1[i] = gcry_int128_to_mpi((__int128_t)pt1.coeffs[i]);
    mpi_pt2[i] = gcry_int128_to_mpi((__int128_t)pt2.coeffs[i]);
  }

  poly_xmy(d, mpi_pt_prod, mpi_pt1, mpi_pt2, NULL, 0);
  for (size_t i=0; i<d; i++)
    pt.coeffs[i] = (double)gcry_mpi_to_int128(mpi_pt_prod[i]);
  //pt.Delta = pt1.Delta*pt2.Delta;
  he_ckks_decode(m, dcd, &pt, &ctx);
  printf("||msg - msg_dec|| = %g\n", blas_dznrmmax_dist(m, msg, dcd));
  free(msg1);
  free(msg2);
#endif
  TEST_END();
}

__attribute__((unused)) static void
test_encrypt_decrypt(int8_t sk[GPQHE_INDCPA_SECRETKEYBYTES],
                     struct ckks_pk *pk)
{
  TEST_BEGIN();
  /* encode */
  double pt[GPQHE_INDCPA_MSGBYTES];
  COMPLEX msg[GPQHE_CKKS_SLOT];
  COMPLEX dec[GPQHE_CKKS_SLOT];
  sample_z01vec(GPQHE_CKKS_SLOT, msg);
  he_ckks_encode(pt, GPQHE_CKKS_SLOT, msg);
  struct ckks_ct ct;
  double diff;

  /* public key encryption */
  TEST_DO("pk enc");
  indcpa_enc_pk(&ct, pt, pk);
#if GPQHE_PMU
  uint64_t t[NTESTS];
  for(uint32_t i=0; i<NTESTS; i++) {
    t[i] = rdtsc();
    indcpa_enc_pk(&ct, pt, pk);
  }
  print_cpuinfo("he_ckks_encrypt_pk", t, NTESTS);
  struct timespec clk[NTESTS];
  for (uint32_t i=0; i<NTESTS; i++){
    clock_gettime(CLOCK_MONOTONIC, &clk[i]);
    indcpa_enc_pk(&ct, pt, pk);
  }
  print_clkinfo("he_ckks_encrypt_pk", clk, NTESTS);
  uint64_t maxrss[NTESTS];
  for(uint32_t i=0; i<NTESTS; i++) {
    maxrss[i] = rdmaxrss();
    indcpa_enc_pk(&ct, pt, pk);
  }
  print_meminfo("he_ckks_encrypt_pk", maxrss, NTESTS);
#endif
  indcpa_dec(pt, &ct, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, dec, pt);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, msg, dec));
  ASSERT(diff<-5);
  printf("log10(||msg'-msg||) = %.15g\n", diff);
  TEST_DONE();

  /* secret key encryption */
  TEST_DO("sk enc");
  indcpa_enc_sk(&ct, pt, sk);
#if GPQHE_PMU
  for(uint32_t i=0; i<NTESTS; i++) {
    t[i] = rdtsc();
    indcpa_enc_sk(&ct, pt, sk);
  }
  print_cpuinfo("he_ckks_encrypt_sk", t, NTESTS);
  for (uint32_t i=0; i<NTESTS; i++){
    clock_gettime(CLOCK_MONOTONIC, &clk[i]);
    indcpa_enc_sk(&ct, pt, sk);
  }
  print_clkinfo("he_ckks_encrypt_sk", clk, NTESTS);
  for(uint32_t i=0; i<NTESTS; i++) {
    maxrss[i] = rdmaxrss();
    indcpa_enc_sk(&ct, pt, sk);
  }
  print_meminfo("he_ckks_encrypt_sk", maxrss, NTESTS);
#endif
  indcpa_dec(pt, &ct, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, dec, pt);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, msg, dec));
  ASSERT(diff<-5);
  printf("log10(||msg'-msg||) = %.15g\n", diff);
  TEST_DONE();
#if GPQHE_PMU
  TEST_DO("decrypt");
  for(uint32_t i=0; i<NTESTS; i++) {
    t[i] = rdtsc();
    indcpa_dec(pt, &ct, sk);
  }
  print_cpuinfo("he_ckks_decrypt", t, NTESTS);
  for (uint32_t i=0; i<NTESTS; i++){
    clock_gettime(CLOCK_MONOTONIC, &clk[i]);
    indcpa_dec(pt, &ct, sk);
  }
  print_clkinfo("he_ckks_decrypt", clk, NTESTS);
  for(uint32_t i=0; i<NTESTS; i++) {
    maxrss[i] = rdmaxrss();
    indcpa_dec(pt, &ct, sk);
  }
  print_meminfo("he_ckks_decrypt", maxrss, NTESTS);
  TEST_DONE();
#endif
  TEST_END();
}

__attribute__((unused)) static void
test_addsub(int8_t sk[GPQHE_INDCPA_SECRETKEYBYTES],
            struct ckks_pk *pk)
{
  TEST_BEGIN();
  /* encode */
  double pt0[GPQHE_INDCPA_MSGBYTES];
  double pt1[GPQHE_INDCPA_MSGBYTES];
  double pt2[GPQHE_INDCPA_MSGBYTES];
  double two[GPQHE_INDCPA_MSGBYTES];
  COMPLEX msg0[GPQHE_CKKS_SLOT];
  COMPLEX msg1[GPQHE_CKKS_SLOT];
  COMPLEX msg2[GPQHE_CKKS_SLOT];
  COMPLEX answer[GPQHE_CKKS_SLOT];
  sample_z01vec(GPQHE_CKKS_SLOT, msg1);
  sample_z01vec(GPQHE_CKKS_SLOT, msg2);
  blas_zxpy(GPQHE_CKKS_SLOT, answer, msg1, msg2);
  he_ckks_encode(pt1, GPQHE_CKKS_SLOT, msg1);
  he_ckks_encode(pt2, GPQHE_CKKS_SLOT, msg2);
  he_ckks_const_pt(two, 2);
  struct ckks_ct ct0, ct1, ct2;
  double diff;

  /* public key encryption */
  TEST_DO("pk encrypted add");
  indcpa_enc_pk(&ct1, pt1, pk);
  indcpa_enc_pk(&ct2, pt2, pk);
#if GPQHE_PMU
  uint64_t t[NTESTS];
  for(uint32_t i=0; i<NTESTS; i++) {
    t[i] = rdtsc();
    he_ckks_add(&ct0, &ct1, &ct2);
  }
  print_cpuinfo("he_ckks_add", t, NTESTS);
  struct timespec clk[NTESTS];
  for (uint32_t i=0; i<NTESTS; i++){
    clock_gettime(CLOCK_MONOTONIC, &clk[i]);
    he_ckks_add(&ct0, &ct1, &ct2);
  }
  print_clkinfo("he_ckks_add", clk, NTESTS);
  uint64_t maxrss[NTESTS];
  for(uint32_t i=0; i<NTESTS; i++) {
    maxrss[i] = rdmaxrss();
    he_ckks_add(&ct0, &ct1, &ct2);
  }
  print_meminfo("he_ckks_add", maxrss, NTESTS);
#endif
  he_ckks_add(&ct0, &ct1, &ct2);
  indcpa_dec(pt0, &ct0, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg0, pt0);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg0));
  ASSERT(diff<-5);
  printf("log10(||msg'-msg||) = %.15g\n", diff);
  TEST_DONE();

  /* secret key encryption */
  TEST_DO("sk encrypted add");
  indcpa_enc_sk(&ct1, pt1, sk);
  indcpa_enc_sk(&ct2, pt2, sk);
#if GPQHE_PMU
  for(uint32_t i=0; i<NTESTS; i++) {
    t[i] = rdtsc();
    he_ckks_add(&ct0, &ct1, &ct2);
  }
  print_cpuinfo("he_ckks_add", t, NTESTS);
  for (uint32_t i=0; i<NTESTS; i++){
    clock_gettime(CLOCK_MONOTONIC, &clk[i]);
    he_ckks_add(&ct0, &ct1, &ct2);
  }
  print_clkinfo("he_ckks_add", clk, NTESTS);
  for(uint32_t i=0; i<NTESTS; i++) {
    maxrss[i] = rdmaxrss();
    he_ckks_add(&ct0, &ct1, &ct2);
  }
  print_meminfo("he_ckks_add", maxrss, NTESTS);
#endif
  he_ckks_add(&ct0, &ct1, &ct2);
  indcpa_dec(pt0, &ct0, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg0, pt0);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg0));
  ASSERT(diff<-5);
  printf("log10(||msg'-msg||) = %.15g\n", diff);
  TEST_DONE();

  /* pk encrypted add plain */
  TEST_DO("pk encrypted add plain");
  indcpa_enc_pk(&ct1, pt1, pk);
#if GPQHE_PMU
  for(uint32_t i=0; i<NTESTS; i++) {
    t[i] = rdtsc();
    he_ckks_add_pt(&ct0, &ct1, pt2);
  }
  print_cpuinfo("he_ckks_add_pt", t, NTESTS);
  for (uint32_t i=0; i<NTESTS; i++){
    clock_gettime(CLOCK_MONOTONIC, &clk[i]);
    he_ckks_add_pt(&ct0, &ct1, pt2);
  }
  print_clkinfo("he_ckks_add_pt", clk, NTESTS);
  for(uint32_t i=0; i<NTESTS; i++) {
    maxrss[i] = rdmaxrss();
    he_ckks_add_pt(&ct0, &ct1, pt2);
  }
  print_meminfo("he_ckks_add_pt", maxrss, NTESTS);
#endif
  he_ckks_add_pt(&ct0, &ct1, pt2);
  indcpa_dec(pt0, &ct0, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg0, pt0);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg0));
  ASSERT(diff<-5);
  printf("log10(||msg'-msg||) = %.15g\n", diff);
  TEST_DONE();

  /* sk encrypted add plain */
  TEST_DO("sk encrypted add plain");
  indcpa_enc_sk(&ct1, pt1, sk);
#if GPQHE_PMU
  for(uint32_t i=0; i<NTESTS; i++) {
    t[i] = rdtsc();
    he_ckks_add(&ct0, &ct1, &ct2);
  }
  print_cpuinfo("he_ckks_add", t, NTESTS);
  for (uint32_t i=0; i<NTESTS; i++){
    clock_gettime(CLOCK_MONOTONIC, &clk[i]);
    he_ckks_add(&ct0, &ct1, &ct2);
  }
  print_clkinfo("he_ckks_add", clk, NTESTS);
  for(uint32_t i=0; i<NTESTS; i++) {
    maxrss[i] = rdmaxrss();
    he_ckks_add(&ct0, &ct1, &ct2);
  }
  print_meminfo("he_ckks_add", maxrss, NTESTS);
#endif
  he_ckks_add_pt(&ct0, &ct1, pt2);
  indcpa_dec(pt0, &ct0, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg0, pt0);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg0));
  ASSERT(diff<-5);
  printf("log10(||msg'-msg||) = %.15g\n", diff);

  /* subtract */
  TEST_DO("pk encrypted sub");
  blas_zxsy(GPQHE_CKKS_SLOT, answer, msg1, msg2);
  indcpa_enc_pk(&ct1, pt1, pk);
  indcpa_enc_pk(&ct2, pt2, pk);
  he_ckks_sub(&ct0, &ct1, &ct2);
  indcpa_dec(pt0, &ct0, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg0, pt0);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg0));
  ASSERT(diff<-5);
  printf("log10(||msg'-msg||) = %.15g\n", diff);
  TEST_DONE();

  /* answer = msg1+2 */
  TEST_DO("msg1+2 (pk enc)");
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    answer[i] = msg1[i] + 2;
  indcpa_enc_pk(&ct1, pt1, pk);
  he_ckks_add_pt(&ct0, &ct1, two);
  indcpa_dec(pt0, &ct0, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg0, pt0);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg0));
  ASSERT(diff<-5);
  printf("log10(||msg'-msg||) = %.15g\n", diff);
  TEST_DONE();

  /* answer = 2-msg2 */
  TEST_DO("2-msg2 (pk enc)");
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    answer[i] = 2-msg2[i];
  indcpa_enc_pk(&ct2, pt2, pk);
  he_ckks_sub_pt(&ct0, &ct2, two);
  he_ckks_neg(&ct0);
  indcpa_dec(pt0, &ct0, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg0, pt0);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg0));
  ASSERT(diff<-5);
  printf("log10(||msg'-msg||) = %.15g\n", diff);
  TEST_DONE();

  TEST_END();
}

__attribute__((unused)) static void
test_multiply(const int8_t sk[GPQHE_INDCPA_SECRETKEYBYTES],
              const struct ckks_pk *pk,
              const struct ckks_pk *rlk)
{
  TEST_BEGIN();
  /* encode */
  double pt0[GPQHE_INDCPA_MSGBYTES];
  double pt1[GPQHE_INDCPA_MSGBYTES];
  double pt2[GPQHE_INDCPA_MSGBYTES];
  COMPLEX msg0[GPQHE_CKKS_SLOT];
  COMPLEX msg1[GPQHE_CKKS_SLOT];
  COMPLEX msg2[GPQHE_CKKS_SLOT];
  COMPLEX answer[GPQHE_CKKS_SLOT];
  sample_z01vec(GPQHE_CKKS_SLOT, msg1);
  sample_z01vec(GPQHE_CKKS_SLOT, msg2);
  blas_zxmy(GPQHE_CKKS_SLOT, answer, msg1, msg2);
  he_ckks_encode(pt1, GPQHE_CKKS_SLOT, msg1);
  he_ckks_encode(pt2, GPQHE_CKKS_SLOT, msg2);
  struct ckks_ct ct0, ct1, ct2;
  double diff;

  /* public key encryption */
  TEST_DO("pk encrypted mult");
  indcpa_enc_pk(&ct1, pt1, pk);
  he_ckks_ct_show_params(&ct1, "ct1 (freshly encrypted)");/* */
  indcpa_enc_pk(&ct2, pt2, pk);
  he_ckks_ct_show_params(&ct1, "ct2 (freshly encrypted)");/* */
#if GPQHE_PMU
  uint64_t t[NTESTS];
  for(uint32_t i=0; i<NTESTS; i++) {
    t[i] = rdtsc();
    he_ckks_mult(&ct0, &ct1, &ct2, rlk);
  }
  print_cpuinfo("he_ckks_mult", t, NTESTS);
  struct timespec clk[NTESTS];
  for (uint32_t i=0; i<NTESTS; i++){
    clock_gettime(CLOCK_MONOTONIC, &clk[i]);
    he_ckks_mult(&ct0, &ct1, &ct2, rlk);
  }
  print_clkinfo("he_ckks_mult", clk, NTESTS);
  uint64_t maxrss[NTESTS];
  for(uint32_t i=0; i<NTESTS; i++) {
    maxrss[i] = rdmaxrss();
    he_ckks_mult(&ct0, &ct1, &ct2, rlk);
  }
  print_meminfo("he_ckks_mult", maxrss, NTESTS);
#endif
  he_ckks_mult(&ct0, &ct1, &ct2, rlk);
  he_ckks_ct_show_params(&ct0, "ct1*ct2");/* */
  he_ckks_rescale(&ct0);
  he_ckks_ct_show_params(&ct0, "ct1*ct2 (rescale)");/* */
  indcpa_dec(pt0, &ct0, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg0, pt0);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg0));
  ASSERT(diff<-5);
  printf("log10(||msg-answer||) = %.15g\n", diff);
  TEST_DONE();

  /* secret key encryption */
  TEST_DO("sk encrypted mult");
  indcpa_enc_sk(&ct1, pt1, sk);
  /* he_ckks_ct_show_params(&ct1, "ct1 (freshly encrypted)"); */
  indcpa_enc_sk(&ct2, pt2, sk);
  /* he_ckks_ct_show_params(&ct1, "ct2 (freshly encrypted)"); */
#if GPQHE_PMU
  for(uint32_t i=0; i<NTESTS; i++) {
    t[i] = rdtsc();
    he_ckks_mult(&ct0, &ct1, &ct2, rlk);
  }
  print_cpuinfo("he_ckks_mult", t, NTESTS);
  for (uint32_t i=0; i<NTESTS; i++){
    clock_gettime(CLOCK_MONOTONIC, &clk[i]);
    he_ckks_mult(&ct0, &ct1, &ct2, rlk);
  }
  print_clkinfo("he_ckks_mult", clk, NTESTS);
  for(uint32_t i=0; i<NTESTS; i++) {
    maxrss[i] = rdmaxrss();
    he_ckks_mult(&ct0, &ct1, &ct2, rlk);
  }
  print_meminfo("he_ckks_mult", maxrss, NTESTS);
#endif
  he_ckks_mult(&ct0, &ct1, &ct2, rlk);
  /* he_ckks_ct_show_params(&ct0, "ct1*ct2"); */
  he_ckks_rescale(&ct0);
  /* he_ckks_ct_show_params(&ct0, "ct1*ct2 (rescale)"); */
  indcpa_dec(pt0, &ct0, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg0, pt0);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg0));
  ASSERT(diff<-5);
  printf("log10(||msg-answer||) = %.15g\n", diff);
  TEST_DONE();

  /* pk encrypted mult plain */
  TEST_DO("pk encrypted plain mult");
  indcpa_enc_pk(&ct1, pt1, pk);
  /* he_ckks_ct_show_params(&ct1, "ct (freshly encrypted)");*/
#if GPQHE_PMU
  for(uint32_t i=0; i<NTESTS; i++) {
    t[i] = rdtsc();
    he_ckks_mult_pt(&ct0, &ct1, pt2);
  }
  print_cpuinfo("he_ckks_mult_pt", t, NTESTS);
  for (uint32_t i=0; i<NTESTS; i++){
    clock_gettime(CLOCK_MONOTONIC, &clk[i]);
    he_ckks_mult_pt(&ct0, &ct1, pt2);
  }
  print_clkinfo("he_ckks_mult_pt", clk, NTESTS);
  for(uint32_t i=0; i<NTESTS; i++) {
    maxrss[i] = rdmaxrss();
    he_ckks_mult_pt(&ct0, &ct1, pt2);
  }
  print_meminfo("he_ckks_mult_pt", maxrss, NTESTS);
#endif
  he_ckks_mult_pt(&ct0, &ct1, pt2);
  he_ckks_rescale(&ct0);
  indcpa_dec(pt0, &ct0, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg0, pt0);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg0));
  ASSERT(diff<-5);
  printf("log10(||msg-answer||) = %.15g\n", diff);
  TEST_DONE();

  /* sk encrypted mult plain */
  TEST_DO("sk encrypted plain mult");
  indcpa_enc_sk(&ct1, pt1, sk);
  /* he_ckks_ct_show_params(&ct1, "ct (freshly encrypted)"); */
  /* he_ckks_ct_show_params(&ct1, "ct (freshly encrypted)");*/
#if GPQHE_PMU
  for(uint32_t i=0; i<NTESTS; i++) {
    t[i] = rdtsc();
    he_ckks_mult_pt(&ct0, &ct1, pt2);
  }
  print_cpuinfo("he_ckks_mult_pt", t, NTESTS);
  for (uint32_t i=0; i<NTESTS; i++){
    clock_gettime(CLOCK_MONOTONIC, &clk[i]);
    he_ckks_mult_pt(&ct0, &ct1, pt2);
  }
  print_clkinfo("he_ckks_mult_pt", clk, NTESTS);
  for(uint32_t i=0; i<NTESTS; i++) {
    maxrss[i] = rdmaxrss();
    he_ckks_mult_pt(&ct0, &ct1, pt2);
  }
  print_meminfo("he_ckks_mult_pt", maxrss, NTESTS);
#endif
  he_ckks_mult_pt(&ct0, &ct1, pt2);
  /* he_ckks_ct_show_params(&ct0, "ct*pt"); */
  he_ckks_rescale(&ct0);
  /* he_ckks_ct_show_params(&ct0, "ct*pt (rescale)"); */
  indcpa_dec(pt0, &ct0, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg0, pt0);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg0));
  ASSERT(diff<-5);
  printf("log10(||msg-answer||) = %.15g\n", diff);
  TEST_DONE();

  /* sk encrypted mult plain (in place) */
  TEST_DO("sk encrypted plain mult (in place)");
  indcpa_enc_sk(&ct0, pt1, sk);
  he_ckks_mult_pt(&ct0, &ct0, pt2);
  he_ckks_rescale(&ct0);
  indcpa_dec(pt0, &ct0, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg0, pt0);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg0));
  ASSERT(diff<-5);
  printf("log10(||msg-answer||) = %.15g\n", diff);
  TEST_DONE();

  TEST_END();
}

__attribute__((unused)) static void
test_conj(const int8_t sk[GPQHE_INDCPA_SECRETKEYBYTES],
          const struct ckks_pk *pk,
          const struct ckks_pk *ck)
{
  TEST_BEGIN();
  /* encode */
  double pt0[GPQHE_INDCPA_MSGBYTES];
  double pt1[GPQHE_INDCPA_MSGBYTES];
  struct ckks_ct ct0, ct1;
  COMPLEX msg0[GPQHE_CKKS_SLOT];
  COMPLEX msg1[GPQHE_CKKS_SLOT];
  COMPLEX answer[GPQHE_CKKS_SLOT];
  sample_z01vec(GPQHE_CKKS_SLOT, msg1);
  he_ckks_encode(pt1, GPQHE_CKKS_SLOT, msg1);
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    answer[i] = conj(msg1[i]);
  double diff;

  /* public key encryption */
  TEST_DO("pk encrypted conj ");
  indcpa_enc_pk(&ct1, pt1, pk);
#if GPQHE_PMU
  uint64_t t[NTESTS];
  for(uint32_t i=0; i<NTESTS; i++) {
    t[i] = rdtsc();
    he_ckks_conj(&ct0, &ct1, ck);
  }
  print_cpuinfo("he_ckks_conj", t, NTESTS);
  struct timespec clk[NTESTS];
  for (uint32_t i=0; i<NTESTS; i++){
    clock_gettime(CLOCK_MONOTONIC, &clk[i]);
    he_ckks_conj(&ct0, &ct1, ck);
  }
  print_clkinfo("he_ckks_conj", clk, NTESTS);
  uint64_t maxrss[NTESTS];
  for(uint32_t i=0; i<NTESTS; i++) {
    maxrss[i] = rdmaxrss();
    he_ckks_conj(&ct0, &ct1, ck);
  }
  print_meminfo("he_ckks_conj", maxrss, NTESTS);
#endif
  he_ckks_conj(&ct0, &ct1, ck);
  indcpa_dec(pt0, &ct0, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg0, pt0);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg0));
  ASSERT(diff<-5);
  printf("log10(||msg'-msg||) = %.15g\n", diff);
  TEST_DONE();

  /* secret key encryption */
  TEST_DO("sk encrypted conj");
  indcpa_enc_sk(&ct1, pt1, sk);
#if GPQHE_PMU
  for(uint32_t i=0; i<NTESTS; i++) {
    t[i] = rdtsc();
    he_ckks_conj(&ct0, &ct1, ck);
  }
  print_cpuinfo("he_ckks_conj", t, NTESTS);
  for (uint32_t i=0; i<NTESTS; i++){
    clock_gettime(CLOCK_MONOTONIC, &clk[i]);
    he_ckks_conj(&ct0, &ct1, ck);
  }
  print_clkinfo("he_ckks_conj", clk, NTESTS);
  for(uint32_t i=0; i<NTESTS; i++) {
    maxrss[i] = rdmaxrss();
    he_ckks_conj(&ct0, &ct1, ck);
  }
  print_meminfo("he_ckks_conj", maxrss, NTESTS);
#endif
  he_ckks_conj(&ct0, &ct1, ck);
  indcpa_dec(pt0, &ct0, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg0, pt0);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg0));
  ASSERT(diff<-5);
  printf("log10(||msg'-msg||) = %.15g\n", diff);
  TEST_DONE();

  TEST_END();
}

__attribute__((unused)) static void
test_rot(const int8_t sk[GPQHE_INDCPA_SECRETKEYBYTES],
         const struct ckks_pk *pk,
         const struct ckks_pk *rk)
{
  TEST_BEGIN();
  /* encode */
  double pt0[GPQHE_INDCPA_MSGBYTES];
  double pt1[GPQHE_INDCPA_MSGBYTES];
  struct ckks_ct ct0, ct1;
  COMPLEX msg0[GPQHE_CKKS_SLOT];
  COMPLEX msg1[GPQHE_CKKS_SLOT];
  COMPLEX answer[GPQHE_CKKS_SLOT];
  sample_z01vec(GPQHE_CKKS_SLOT, msg1);
  double diff;

  TEST_DO("pk encrypted rot all");
  for (size_t rot=0; rot<GPQHE_CKKS_SLOT; rot++)
  {
    blas_zrot(GPQHE_CKKS_SLOT, answer, rot, msg1);
    he_ckks_encode(pt1, GPQHE_CKKS_SLOT, msg1);
    indcpa_enc_pk(&ct1, pt1, pk);
    he_ckks_rot(&ct0, &ct1, rot, &rk[rot]);
    indcpa_dec(pt0, &ct0, sk);
    he_ckks_decode(GPQHE_CKKS_SLOT, msg0, pt0);
    diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg0));
    ASSERT(diff<-5);
    printf("rot=%ld: log10(||msg'-msg||) = %.15g\n", rot, diff);
  }
  TEST_DONE();

  /* pk encrypted */
  TEST_DO("pk encrypted rot");
  blas_zrot(GPQHE_CKKS_SLOT, answer, 1, msg1);
  he_ckks_encode(pt1, GPQHE_CKKS_SLOT, msg1);
  indcpa_enc_pk(&ct1, pt1, pk);
#if GPQHE_PMU
  uint64_t t[NTESTS];
  for(uint32_t i=0; i<NTESTS; i++) {
    t[i] = rdtsc();
    he_ckks_rot(&ct0, &ct1, 1, &rk[1]);
  }
  print_cpuinfo("he_ckks_rot", t, NTESTS);
  struct timespec clk[NTESTS];
  for (uint32_t i=0; i<NTESTS; i++){
    clock_gettime(CLOCK_MONOTONIC, &clk[i]);
    he_ckks_rot(&ct0, &ct1, 1, &rk[1]);
  }
  print_clkinfo("he_ckks_rot", clk, NTESTS);
  uint64_t maxrss[NTESTS];
  for(uint32_t i=0; i<NTESTS; i++) {
    maxrss[i] = rdmaxrss();
    he_ckks_rot(&ct0, &ct1, 1, &rk[1]);
  }
  print_meminfo("he_ckks_rot", maxrss, NTESTS);
#endif
  he_ckks_rot(&ct0, &ct1, 1, &rk[1]);
  indcpa_dec(pt0, &ct0, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg0, pt0);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg0));
  ASSERT(diff<-5);
  printf("rot=%ld: log10(||msg'-msg||) = %.15g\n", 1l, diff);
  TEST_DONE();

  /* sk encrypted */
  TEST_DO("sk encrypted");
  indcpa_enc_sk(&ct1, pt1, sk);
#if GPQHE_PMU
  for(uint32_t i=0; i<NTESTS; i++) {
    t[i] = rdtsc();
    he_ckks_rot(&ct0, &ct1, 1, &rk[1]);
  }
  print_cpuinfo("he_ckks_rot", t, NTESTS);
  for (uint32_t i=0; i<NTESTS; i++){
    clock_gettime(CLOCK_MONOTONIC, &clk[i]);
    he_ckks_rot(&ct0, &ct1, 1, &rk[1]);
  }
  print_clkinfo("he_ckks_rot", clk, NTESTS);
  for(uint32_t i=0; i<NTESTS; i++) {
    maxrss[i] = rdmaxrss();
    he_ckks_rot(&ct0, &ct1, 1, &rk[1]);
  }
  print_meminfo("he_ckks_rot", maxrss, NTESTS);
#endif
  he_ckks_rot(&ct0, &ct1, 1, &rk[1]);
  indcpa_dec(pt0, &ct0, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg0, pt0);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg0));
  ASSERT(diff<-5);
  printf("rot=%ld: log10(||msg'-msg||) = %.15g\n", 1l, diff);
  TEST_DONE();

  TEST_END();
}

__attribute__((unused)) static void
test_gemv(const int8_t sk[GPQHE_INDCPA_SECRETKEYBYTES],
          const struct ckks_pk *pk,
          const struct ckks_pk *rk)
{
  TEST_BEGIN();
  /* encode */
  double pt_Av[GPQHE_INDCPA_MSGBYTES];
  double pt_v [GPQHE_INDCPA_MSGBYTES];
  struct ckks_ct ct_Av, ct_v;
  COMPLEX Av[GPQHE_CKKS_SLOT];
  COMPLEX A[(GPQHE_CKKS_SLOT)*(GPQHE_CKKS_SLOT)];
  COMPLEX v[GPQHE_CKKS_SLOT];
  COMPLEX answer[GPQHE_CKKS_SLOT];
  sample_z01vec(GPQHE_CKKS_SLOT*GPQHE_CKKS_SLOT, A);
  sample_z01vec(GPQHE_CKKS_SLOT, v);
  double diff;
  blas_zgemv(GPQHE_CKKS_SLOT, GPQHE_CKKS_SLOT, answer, A, v);
  he_ckks_encode(pt_v, GPQHE_CKKS_SLOT, v);

  /* pk encrypted */
  TEST_DO("pk encrypted");
  indcpa_enc_pk(&ct_v, pt_v, pk);
#if GPQHE_PMU
  uint64_t t[NTESTS];
  for (uint32_t i=0; i<NTESTS; i++){
    t[i] = rdtsc();
    he_ckks_gemv(&ct_Av, A, &ct_v, rk);
  }
  print_cpuinfo("he_ckks_gemv", t, NTESTS);
  struct timespec clk[NTESTS];
  for (uint32_t i=0; i<NTESTS; i++){
    clock_gettime(CLOCK_MONOTONIC, &clk[i]);
    he_ckks_gemv(&ct_Av, A, &ct_v, rk);
  }
  print_clkinfo("he_ckks_gemv", clk, NTESTS);
  uint64_t maxrss[NTESTS];
  for (uint32_t i=0; i<NTESTS; i++){
    maxrss[i] = rdmaxrss();
    he_ckks_gemv(&ct_Av, A, &ct_v, rk);
  }
  print_meminfo("he_ckks_gemv", maxrss, NTESTS);
#endif
  he_ckks_gemv(&ct_Av, A, &ct_v, rk);
  indcpa_dec(pt_Av, &ct_Av, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, Av, pt_Av);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, Av));
  ASSERT(diff<-5);
  printf("log10(||msg'-msg||) = %.15g\n", diff);
  TEST_DONE();

  /* sk encrypted */
  TEST_DO("sk encrypted");
  indcpa_enc_sk(&ct_v, pt_v, sk);
#if GPQHE_PMU
  for (uint32_t i=0; i<NTESTS; i++){
    t[i] = rdtsc();
    he_ckks_gemv(&ct_Av, A, &ct_v, rk);
  }
  print_cpuinfo("he_ckks_gemv", t, NTESTS);
  for (uint32_t i=0; i<NTESTS; i++){
    clock_gettime(CLOCK_MONOTONIC, &clk[i]);
    he_ckks_gemv(&ct_Av, A, &ct_v, rk);
  }
  print_clkinfo("he_ckks_gemv", clk, NTESTS);
  for (uint32_t i=0; i<NTESTS; i++){
    maxrss[i] = rdmaxrss();
    he_ckks_gemv(&ct_Av, A, &ct_v, rk);
  }
  print_meminfo("he_ckks_gemv", maxrss, NTESTS);
#endif
  he_ckks_gemv(&ct_Av, A, &ct_v, rk);
  indcpa_dec(pt_Av, &ct_Av, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, Av, pt_Av);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, Av));
  ASSERT(diff<-5);
  printf("log10(||msg'-msg||) = %.15g\n", diff);
  TEST_DONE();

  TEST_END();
}

/*
 * +-----------+       +-----------+
 * |ct, ct_conj|======>| ct0 | ct1 |
 * +-----------+       +-----------+
 *      ^ |                |    |
 *      | V                V    V
 * +-----------+         +--+  +--+
 * |     pt    |         |pt|  |pt|
 * +-----------+         +--+  +--+
 *      ^ |                |    |
 *      | |                V    V
 *      | |            +-------------+
 *      | V            | dcdr | dcdi |
 * +-----------+       +-------------+
 * |    msg    |       | msgr | msgi |
 * +-----------+       +-------------+
 */
__attribute__((unused)) static void
test_coeff_to_slot(const int8_t sk[GPQHE_INDCPA_SECRETKEYBYTES],
                   const struct ckks_pk *pk,
                   const struct ckks_pk *ck,
                   const struct ckks_pk *rk)
{
  TEST_BEGIN();

  double pt[GPQHE_INDCPA_MSGBYTES];
  COMPLEX msg [GPQHE_CKKS_SLOT];
  COMPLEX msgr[GPQHE_CKKS_SLOT]; /* msg real part */
  COMPLEX msgi[GPQHE_CKKS_SLOT]; /* msg imag part */
  COMPLEX dcdr[GPQHE_CKKS_SLOT];
  COMPLEX dcdi[GPQHE_CKKS_SLOT];
  struct ckks_ct ct, ct0, ct1;
  double diff;
  sample_z01vec(GPQHE_CKKS_SLOT, msg);
  he_ckks_encode(pt, GPQHE_CKKS_SLOT, msg);
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
  {
    msgr[i] = pt[i]/GPQHE_P;
    msgi[i] = pt[i+GPQHE_CKKS_SLOT]/GPQHE_P;
  }
  indcpa_enc_pk(&ct, pt, pk);
  he_ckks_coeff_to_slot(&ct0, &ct1, &ct, ck, rk);
  indcpa_dec(pt, &ct0, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, dcdr, pt);
  indcpa_dec(pt, &ct1, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, dcdi, pt);

  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, dcdr, msgr));
  ASSERT(diff<-5);
  printf("log10(||dcdr-msgr||) = %.20g\n", diff);

  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, dcdi, msgi));
  ASSERT(diff<-5);
  printf("log10(||dcdi-msgi||) = %.20g\n", diff);

  TEST_END();
}

__attribute__((unused)) static void
test_slot_to_coeff(const int8_t sk[GPQHE_INDCPA_SECRETKEYBYTES],
                   const struct ckks_pk *pk,
                   const struct ckks_pk *ck,
                   const struct ckks_pk *rk)
{
  TEST_BEGIN();
  double pt [GPQHE_INDCPA_MSGBYTES];
  double pt0[GPQHE_INDCPA_MSGBYTES];
  double pt1[GPQHE_INDCPA_MSGBYTES];
  double ptn[GPQHE_INDCPA_MSGBYTES];
  double pt_ans0[GPQHE_INDCPA_MSGBYTES];
  double pt_ans1[GPQHE_INDCPA_MSGBYTES];
  struct ckks_ct ct, ct0, ct1;
  COMPLEX msg [GPQHE_CKKS_SLOT];
  COMPLEX msg0[GPQHE_CKKS_SLOT];
  COMPLEX msg1[GPQHE_CKKS_SLOT];
  COMPLEX msg_ans0[GPQHE_CKKS_SLOT];
  COMPLEX msg_ans1[GPQHE_CKKS_SLOT];
  double diff;

  sample_z01vec(GPQHE_CKKS_SLOT, msg);
  he_ckks_encode(pt, GPQHE_CKKS_SLOT, msg);
  indcpa_enc_pk(&ct, pt, pk);

  he_ckks_coeff_to_slot(&ct0, &ct1, &ct, ck, rk);
  he_ckks_slot_to_coeff(&ct, &ct0, &ct1, rk);
  indcpa_dec(ptn, &ct, sk);
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
  {
    pt_ans0[i] = ptn[i]/GPQHE_P;
    pt_ans1[i] = ptn[i+GPQHE_CKKS_SLOT]/GPQHE_P;
  }
  printf("ptn\n");
  for (size_t i=0; i<GPQHE_INDCPA_MSGBYTES; i++)
    printf("% f\n", ptn[i]);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg_ans0, pt_ans0);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg_ans1, pt_ans1);

  indcpa_dec(pt0, &ct0, sk);
  indcpa_dec(pt1, &ct1, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg0, pt0);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg1, pt1);
  printf("msg0                  msg_ans0\n");
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    printf("% f % fj  % f % fj\n",
      creal(msg0[i]), cimag(msg0[i]),
      creal(msg_ans0[i]), cimag(msg_ans0[i]));
  printf("msg1                  msg_ans1\n");
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    printf("% f % fj  % f % fj\n",
      creal(msg1[i]), cimag(msg1[i]),
      creal(msg_ans1[i]), cimag(msg_ans1[i]));

  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, msg_ans0, msg0));
  ASSERT(diff<-5);
  printf("log10(||msg_ans0-msg0||) = %.20g\n",diff);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, msg_ans1, msg1));
  ASSERT(diff<-5);
  printf("log10(||msg_ans1-msg1||) = %.20g\n",diff);
  TEST_END();
}

__attribute__((unused)) static void
test_exp(const int8_t sk[GPQHE_INDCPA_SECRETKEYBYTES],
         const struct ckks_pk *pk __attribute__((unused)),
         const struct ckks_pk *rlk)
{
  TEST_BEGIN();
  COMPLEX msg[GPQHE_CKKS_SLOT];
  COMPLEX msg_exp[GPQHE_CKKS_SLOT];
  COMPLEX answer[GPQHE_CKKS_SLOT];
  COMPLEX const_pt[GPQHE_CKKS_SLOT];
  double pt[GPQHE_INDCPA_MSGBYTES];
  struct ckks_ct ct, a_ct, ct_exp;
  COMPLEX a;
  size_t iteration = 5;
  double diff;
  sample_z01vec(GPQHE_CKKS_SLOT, msg);
  /* encrypt to ct */
  he_ckks_encode(pt, GPQHE_CKKS_SLOT, msg);
  indcpa_enc_sk(&ct, pt, sk);
  he_ckks_ct_show_params(&ct, "ct (fresh)");

  /* a=1/Delta */
  TEST_DO("a=1/Delta, msg=exp(a*ct)");
  a = 1.0/GPQHE_P;
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    answer[i] = cexp(a*msg[i]);
  he_ckks_const_pt(pt, a);
  he_ckks_mult_pt(&a_ct, &ct, pt);
  he_ckks_rescale(&a_ct); /* l-1 */
  he_ckks_exp(&ct_exp, &a_ct, rlk);
#if GPQHE_PMU
  uint64_t cpu[NTESTS];
  for (uint32_t i=0; i<NTESTS; i++){
    cpu[i] = rdtsc();
    he_ckks_exp(&ct_exp, &a_ct, rlk);
  }
  print_cpuinfo("he_ckks_exp", cpu, NTESTS);
  struct timespec clk[NTESTS];
  for (uint32_t i=0; i<NTESTS; i++){
    clock_gettime(CLOCK_MONOTONIC, &clk[i]);
    he_ckks_exp(&ct_exp, &a_ct, rlk);
  }
  print_clkinfo("he_ckks_exp", clk, NTESTS);
  uint64_t maxrss[NTESTS];
  for (uint32_t i=0; i<NTESTS; i++){
    maxrss[i] = rdmaxrss();
    he_ckks_exp(&ct_exp, &a_ct, rlk);
  }
  print_meminfo("he_ckks_exp", maxrss, NTESTS);
#endif
  he_ckks_ct_show_params(&ct_exp, "exp(a*ct)");
  indcpa_dec(pt, &ct_exp, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg_exp, pt);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg_exp));
  ASSERT(diff<-5);
  printf("log10(||msg'-msg||) = %.20g\n", diff);
  TEST_DONE();

  /* a=2Pi*i/Delta */
  TEST_DO("a=2Pi*i/Delta, msg=exp(a*ct)");
  a = 2*M_PI*I/GPQHE_P;
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    answer[i] = cexp(a*msg[i]);
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    const_pt[i] = a;///(1<<ITERATION);
  //he_ckks_const_pt(pt, const_pt);
  he_ckks_encode(pt, GPQHE_CKKS_SLOT, const_pt);
  he_ckks_mult_pt(&a_ct, &ct, pt);
  he_ckks_rescale(&a_ct); /* l-1 */
  he_ckks_exp(&ct_exp, &a_ct, rlk);
  he_ckks_ct_show_params(&ct_exp, "exp(a*ct)");
  indcpa_dec(pt, &ct_exp, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg_exp, pt);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg_exp));
  ASSERT(diff<-5);
  printf("log10(||msg'-msg||) = %.20g\n", diff);
  TEST_DONE();

  /* a=2Pi*i/(2^r) */
  TEST_DO("a=2Pi*i/(2^r), msg=exp(a*ct)^(2^r) (r=5)");
  a = 2*M_PI*I;
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    answer[i] = cexp(a*msg[i]);
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    const_pt[i] = a/(1<<iteration);
  he_ckks_encode(pt, GPQHE_CKKS_SLOT, const_pt);
  he_ckks_mult_pt(&a_ct, &ct, pt);
  he_ckks_rescale(&a_ct); /* l-1 */
  he_ckks_exp(&ct_exp, &a_ct, rlk);
  for (size_t i=0; i<iteration; i++){
    he_ckks_mult(&ct_exp, &ct_exp, &ct_exp, rlk);
    he_ckks_rescale(&ct_exp);
  }
  he_ckks_ct_show_params(&ct_exp, "exp(a*ct)");
  indcpa_dec(pt, &ct_exp, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg_exp, pt);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg_exp));
  ASSERT(diff<-5);
  printf("log10(||msg'-msg||) = %.20g\n", diff);
  TEST_DONE();

  /* a=2*Pi/(2^r) */
  TEST_DO("a=2*Pi/(2^r), msg=exp(a*ct)^(2^r) (r=5)");
  a = 2*M_PI;
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    answer[i] = cexp(a*msg[i]);
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    const_pt[i] = a/(1<<iteration);
  he_ckks_encode(pt, GPQHE_CKKS_SLOT, const_pt);
  he_ckks_mult_pt(&a_ct, &ct, pt);
  he_ckks_rescale(&a_ct); /* l-1 */
  he_ckks_exp(&ct_exp, &a_ct, rlk);
  for (size_t i=0; i<iteration; i++){
    he_ckks_mult(&ct_exp, &ct_exp, &ct_exp, rlk);
    he_ckks_rescale(&ct_exp);
  }
  he_ckks_ct_show_params(&ct_exp, "exp(a*ct)");
  indcpa_dec(pt, &ct_exp, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg_exp, pt);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg_exp));
  ASSERT(diff<-5);
  printf("log10(||msg'-msg||) = %.20g\n", diff);
  TEST_DONE();

  /* a=2Pi*i/(Delta*2^r) */
  TEST_DO("a=2Pi*i/(Delta*2^r), msg=exp(a*ct)^(2^r) (r=5)");
  a = 2*M_PI*I/GPQHE_P;
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    answer[i] = cexp(a*msg[i]);
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    const_pt[i] = a/(1<<iteration);
  he_ckks_encode(pt, GPQHE_CKKS_SLOT, const_pt);
  he_ckks_mult_pt(&a_ct, &ct, pt);
  he_ckks_rescale(&a_ct); /* l-1 */
  he_ckks_exp(&ct_exp, &a_ct, rlk);
  for (size_t i=0; i<iteration; i++){
    he_ckks_mult(&ct_exp, &ct_exp, &ct_exp, rlk);
    he_ckks_rescale(&ct_exp);
  }
  he_ckks_ct_show_params(&ct_exp, "exp(a*ct)");
  indcpa_dec(pt, &ct_exp, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg_exp, pt);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg_exp));
  ASSERT(diff<-5);
  printf("log10(||msg'-msg||) = %.20g\n", diff);
  TEST_DONE();

  /* a=2Pi/(Delta*2^r) */
  TEST_DO("a=2Pi/(Delta*2^r), msg=exp(a*ct)^(2^r) (r=5)");
  a = 2*M_PI/GPQHE_P;
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    answer[i] = cexp(a*msg[i]);
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    const_pt[i] = a/(1<<iteration);
  he_ckks_encode(pt, GPQHE_CKKS_SLOT, const_pt);
  he_ckks_mult_pt(&a_ct, &ct, pt);
  he_ckks_rescale(&a_ct); /* l-1 */
  he_ckks_exp(&ct_exp, &a_ct, rlk);
  for (size_t i=0; i<iteration; i++){
    he_ckks_mult(&ct_exp, &ct_exp, &ct_exp, rlk);
    he_ckks_rescale(&ct_exp);
  }
  he_ckks_ct_show_params(&ct_exp, "exp(a*ct)");
  indcpa_dec(pt, &ct_exp, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg_exp, pt);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg_exp));
  ASSERT(diff<-5);
  printf("log10(||msg'-msg||) = %.20g\n", diff);
  TEST_DONE();

  TEST_END();
}

__attribute__ ((unused)) static void
test_rlsin(const int8_t sk[GPQHE_INDCPA_SECRETKEYBYTES],
           const struct ckks_pk *pk __attribute__((unused)),
           const struct ckks_pk *rlk,
           const struct ckks_pk *ck)
{
  TEST_BEGIN();
  COMPLEX msg[GPQHE_CKKS_SLOT];
  COMPLEX msg_rlsin[GPQHE_CKKS_SLOT];
  COMPLEX answer[GPQHE_CKKS_SLOT];
  double pt[GPQHE_INDCPA_MSGBYTES];
  struct ckks_ct ct, ct_rlsin;
  double a;
  double diff;
  sample_z01vec(GPQHE_CKKS_SLOT, msg);
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    answer[i] = msg[i]/GPQHE_P;
  /* encrypt to ct */
  he_ckks_encode(pt, GPQHE_CKKS_SLOT, msg);
  indcpa_enc_sk(&ct, pt, sk);

  TEST_DO("a=2*Pi");
  a=2*M_PI;
  he_ckks_rlsin(&ct_rlsin, &ct, a, rlk, ck);
  indcpa_dec(pt, &ct_rlsin, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg_rlsin, pt);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg_rlsin));
  ASSERT(diff<-5);
  printf("log10(||msg'-msg||) = %.20g\n", diff);
  TEST_DONE();

  TEST_DO("a=2*Pi/GPQHE_P");
  a=2*M_PI/GPQHE_P;
  he_ckks_rlsin(&ct_rlsin, &ct, a, rlk, ck);
  indcpa_dec(pt, &ct_rlsin, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg_rlsin, pt);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg_rlsin));
  ASSERT(diff<-5);
  printf("log10(||msg'-msg||) = %.20g\n", diff);
  TEST_DONE();

  TEST_END();
}

#if 0
__attribute__((unused)) static void
test_he_ckks_sin(ckks_sk *sk, ckks_pk *pk, ckks_evk *evk, ckks_ctx *ctx)
{
  size_t msg_len = GPQHE_N/2;
  COMPLEX *msg     = (COMPLEX *)calloc(msg_len, sizeof(COMPLEX));
  COMPLEX *msg_sin = (COMPLEX *)calloc(msg_len, sizeof(COMPLEX));
  COMPLEX *answer  = (COMPLEX *)calloc(msg_len, sizeof(COMPLEX));
  ckks_pt pt;
  ckks_ct ct, ct_sin;
  double a = 2*M_PI/GPQHE_P;
  sample_z01vec(msg_len, msg);
  /* encrypt to ct */
  he_ckks_encode(&pt, msg_len, msg, ctx);
  he_ckks_encrypt_pk(&ct, &pt, pk, ctx);
  he_ckks_ct_show_params(&ct, "ct (fresh)");
  /* calculate sin(a*ct) */
  for (size_t i=0; i<msg_len; i++)
    answer[i] = csin(a*msg[i]);
  he_ckks_encode(&pt, msg_len, msg, ctx);
  he_ckks_encrypt_pk(&ct, &pt, pk, ctx);
  he_ckks_const_pt(&pt, a);
  he_ckks_mult_pt(&ct, &ct, &pt, ctx);
  he_ckks_rescale(&ct, ctx);
  he_ckks_sin(&ct_sin, &ct, &evk->relin_key, ctx);
  he_ckks_ct_show_params(&ct_sin, "sin(a*ct)");
  he_ckks_decrypt(&pt, &ct_sin, sk, ctx);
  he_ckks_decode(msg_len, msg_sin, &pt, ctx);
  printf("||sin(a*ct)-answer|| = %.20g\n", blas_dznrmmax_dist(msg_len, answer, msg_sin));
  printf("\033[1m%s pass.\033[0m\n", __func__);
}
#endif

__attribute__((unused)) static void
test_bootstrap(const int8_t sk[GPQHE_INDCPA_SECRETKEYBYTES],
               const struct ckks_pk *pk __attribute__((unused)),
               const struct ckks_pk *rlk,
               const struct ckks_pk *ck,
               const struct ckks_pk *rk)
{
  TEST_BEGIN();

  COMPLEX msg[GPQHE_CKKS_SLOT];
  COMPLEX msg_bootstrap[GPQHE_CKKS_SLOT];
  double pt[GPQHE_INDCPA_MSGBYTES];
  struct ckks_ct ct;
  double diff;
  sample_z01vec(GPQHE_CKKS_SLOT, msg);
  /* encrypt to ct */
  he_ckks_encode(pt, GPQHE_CKKS_SLOT, msg);
  indcpa_enc_sk(&ct, pt, sk);

  he_ckks_ct_show_params(&ct, "freshly encrypted");
  while (ct.l>2)
    he_ckks_moddown(&ct);
  he_ckks_ct_show_params(&ct, "by now should do bootstrap");

  TEST_DO("bootstrap");
  he_ckks_bootstrap(&ct, rlk, ck, rk);
  indcpa_dec(pt, &ct, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg_bootstrap, pt);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, msg, msg_bootstrap));
  ASSERT(diff<-5);
  printf("log10(||msg'-msg||) = %.20g\n", diff);
  he_ckks_ct_show_params(&ct, "after bootstrap");
  TEST_DONE();

  TEST_END();
}

__attribute__((unused)) static void
test_inv(const int8_t sk[GPQHE_INDCPA_SECRETKEYBYTES],
         const struct ckks_pk *pk __attribute__((unused)),
         const struct ckks_pk *rlk)
{
  TEST_BEGIN();

  double pt[GPQHE_INDCPA_MSGBYTES];
  struct ckks_ct ct, ct_inv;
  COMPLEX msg[GPQHE_CKKS_SLOT];
  COMPLEX answer[GPQHE_CKKS_SLOT];
  double diff;
  sample_z01vec(GPQHE_CKKS_SLOT, msg);
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
  {
    msg[i] = creal(msg[i])==0? 0.5 : creal(msg[i]);
    answer[i] = 1/msg[i];
  }
  TEST_DO("inv");
  he_ckks_encode(pt, GPQHE_CKKS_SLOT, msg);
  indcpa_enc_pk(&ct, pt, pk);
  he_ckks_inv(&ct_inv, &ct, rlk, 7);
  he_ckks_ct_show_params(&ct_inv, "inv(ct)");
  indcpa_dec(pt, &ct_inv, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg, pt);
  diff = log10(blas_dznrmmax_dist(GPQHE_CKKS_SLOT, answer, msg));
  ASSERT(diff<-5);
  printf("pk encrypted inv: ||msg-answer|| = %.15g\n", diff);
  TEST_DONE();

  TEST_END();
}

__attribute__((unused)) static void
test_sqrt(const int8_t sk[GPQHE_INDCPA_SECRETKEYBYTES],
          const struct ckks_pk *pk __attribute__((unused)),
          const struct ckks_pk *rlk)
{
  TEST_BEGIN();

  size_t iter=1;
  double pt[GPQHE_INDCPA_MSGBYTES];
  struct ckks_ct ct, ct_sqrt;
  COMPLEX msg[GPQHE_CKKS_SLOT];
  COMPLEX msg_sqrt[GPQHE_CKKS_SLOT];
  COMPLEX an[GPQHE_CKKS_SLOT];
  COMPLEX bn[GPQHE_CKKS_SLOT];
  COMPLEX an_ret[GPQHE_CKKS_SLOT];
  COMPLEX bn_ret[GPQHE_CKKS_SLOT];
  COMPLEX tmp_ret[GPQHE_CKKS_SLOT];
  double diff;

  TEST_DO("compare an/bn algorithm with sqrt from libc");
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
  {
    msg[i]=(double)rand()/RAND_MAX;
    msg_sqrt[i] = sqrt(msg[i]);
    an[i] = msg[i];
    bn[i] = msg[i]-1;
    for (size_t n=0; n<iter; n++)
    {
      an[i] = an[i]*(1-bn[i]/2);
      //bn[i] = bn[i]*bn[i]*(bn[i]-3)/4;
    }
  }
  printf("msg[0]=%f+%fj\n", creal(msg[0]), cimag(msg[0]));
  printf("an[0]=%f+%fj\n", creal(an[0]), cimag(an[0]));
  printf("1-bn[0]/2=%f+%fj\n", creal(1-bn[0]/2), cimag(1-bn[0]/2));
#if 0
  diff = blas_dznrmmax_dist(GPQHE_CKKS_SLOT, msg, an);
  ASSERT(diff<1e-5);
  printf("plain: ||msg-answer|| = %.15g\n", diff);
#endif
  TEST_DONE();

  TEST_DO("test HE sqrt");
  he_ckks_encode(pt, GPQHE_CKKS_SLOT, msg);
  indcpa_enc_pk(&ct, pt, pk);
#if 0
  he_ckks_sqrt(&ct_sqrt, &ct, rlk);
#else
  COMPLEX num   [GPQHE_CKKS_SLOT];
  double one    [GPQHE_INDCPA_MSGBYTES];
  double three  [GPQHE_INDCPA_MSGBYTES];
  double half   [GPQHE_INDCPA_MSGBYTES];
  double quarter[GPQHE_INDCPA_MSGBYTES];
  struct ckks_ct tmp, an_, bn_;
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    num[i]=1;
  he_ckks_encode(one, GPQHE_CKKS_SLOT, num);
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    num[i]=3;
  he_ckks_encode(three, GPQHE_CKKS_SLOT, num);
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    num[i]=0.5;
  he_ckks_encode(half, GPQHE_CKKS_SLOT, num);
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    num[i]=0.25;
  he_ckks_encode(quarter, GPQHE_CKKS_SLOT, num);
  he_ckks_ct_copy(&an_, &ct);           /* an = ct */
  he_ckks_sub_pt(&bn_, &ct, one); /* bn = ct-1 */
  for (size_t n=0; n<iter; n++)
  {
    /* an */
    he_ckks_mult_pt(&tmp, &bn_, half); /* bn/2 */
    he_ckks_rescale(&tmp);
    he_ckks_sub_pt(&tmp, &tmp, one);  /* bn/2-1 */
    he_ckks_neg(&tmp);                 /* tmp = 1-bn/2 */
    he_ckks_moddown(&an_);
    he_ckks_mult(&an_, &an_, &an_, rlk); /* an = an*(1-bn/2) */
    he_ckks_rescale(&an_);
#if 0
    /* bn */
    he_ckks_sub_pt(&tmp, &bn_, three);
    he_ckks_mult_pt(&tmp, &tmp, quarter);   /* tmp = (bn-3)/4 */
    he_ckks_rescale(&tmp);
    he_ckks_mult(&bn_, &bn_, &bn_, rlk);  /* bn = bn*bn */
    he_ckks_rescale(&bn_);
    he_ckks_mult(&bn_, &bn_, &tmp, rlk); /* bn = (bn*bn)*((bn-3)/4) */
    he_ckks_rescale(&bn_);
#endif
  }
  he_ckks_ct_copy(&ct_sqrt, &an_);
#endif
  he_ckks_ct_show_params(&ct_sqrt, "sqrt(ct)");
  puts("");

  indcpa_dec(pt, &ct_sqrt, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, msg_sqrt, pt);
  printf("an[0]=%f+%fj\n", creal(an[0]), cimag(an[0]));

  indcpa_dec(pt, &an_, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, an_ret, pt);
  printf("an_ret[0]=%f+%fj\n", creal(an_ret[0]), cimag(an_ret[0]));

  indcpa_dec(pt, &bn_, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, bn_ret, pt);

  indcpa_dec(pt, &tmp, sk);
  he_ckks_decode(GPQHE_CKKS_SLOT, tmp_ret, pt);
  printf("tmp_ret[0]=%f+%fj\n", creal(tmp_ret[0]), cimag(tmp_ret[0]));

  diff = blas_dznrmmax_dist(GPQHE_CKKS_SLOT, an_ret, an);
  printf("||an_ret-an||=%.15g\n",diff);
  diff = blas_dznrmmax_dist(GPQHE_CKKS_SLOT, bn_ret, bn);
  printf("||bn_ret-bn||=%.15g\n",diff);

  diff = blas_dznrmmax_dist(GPQHE_CKKS_SLOT, msg_sqrt, an);
  ASSERT(diff<1e-5);
  printf("pk encrypted sqrt: ||msg-answer|| = %.15g\n", diff);
  TEST_DONE();

  TEST_END();
}

#if 0
__attribute__((unused)) static void
actuate(const char *name, double *x0, double *uopt,
  struct controller_ctx *controller)
{
  uint32_t n = GPQHE_DIMX;
  uint32_t m = GPQHE_DIMU;
  uint32_t N = GPQHE_HORIZON;
  int openfd;
  ssize_t readfd;
  double x[GPQHE_DIMX] = {0};
  double u[GPQHE_DIMU] = {0};
  /* load matrices */
  double A[GPQHE_DIMX*GPQHE_DIMX] = {0};
  openfd = open("precomp/A", O_RDONLY, 0000);
  readfd = read(openfd, A, n*n*sizeof(double));
  double B[GPQHE_DIMX*GPQHE_DIMU] = {0};
  openfd = open("precomp/B", O_RDONLY, 0000);
  readfd = read(openfd, B, n*m*sizeof(double));
  close(readfd);
  close(openfd);
  /* state evolution */
  double record[(GPQHE_DIMX+GPQHE_DIMU)*(GPQHE_HORIZON+1)] = {0};
  x[0]=x0[0],x[1]=x0[1];
  for (uint32_t i=0; i<N; i++)
  {
    for (uint32_t j=0; j<n; j++)
      record[i*(n+m)+j]=x[j];
    for (uint32_t j=0; j<m; j++)
      record[i*(n+m)+n+j]=u[j]=uopt[i*m+j];
    double Ax[GPQHE_DIMX] = {0};
    cblas_dgemv(CblasRowMajor,CblasNoTrans,n,n,1,A,n,x,1,1,Ax,1);
    double Bu[GPQHE_DIMX] = {0};
    cblas_dgemv(CblasRowMajor,CblasNoTrans,n,m,1,B,m,u,1,1,Bu,1);
    cblas_dcopy(n,Ax,1,x,1);
    cblas_daxpy(n,1,Bu,1,x,1);
  }
  for (uint32_t j=0; j<n; j++)
    record[N*(n+m)+j]=x[j];
  openfd = open(name, O_WRONLY | O_CREAT, 0664);
  ssize_t writefd = write(openfd, record, (n+m)*(N+1)*sizeof(double));
  close(writefd);
  close(openfd);
}

__attribute__((unused)) static void
test_finite_horizon_optimal_input(int argc, const char *argv[],
  ckks_sk *sk, ckks_pk *pk, ckks_evk *evk, ckks_ctx *ctx,
  struct controller_ctx *controller)
{
  TEST_BEGIN();
  uint32_t n = GPQHE_DIMX;
  uint32_t m = GPQHE_DIMU;
  uint32_t N = GPQHE_HORIZON;
  char name[255];
  /* initial state */
  argc--;
  argv++;
  assert(argc==GPQHE_DIMX);
  double x0[GPQHE_DIMX];
  for (size_t i=0; i<GPQHE_DIMX; i++)
    x0[i]=atof(argv[i]);
  /* uopt */
  double uopt[GPQHE_HORIZON*GPQHE_DIMU] = {0};
#if GPQHE_PMU==1
  uint64_t t[NTESTS];
  for (uint32_t i=0; i<NTESTS; i++){
    t[i] = rdtsc();
    memset(uopt,0x00,GPQHE_HORIZON*GPQHE_DIMU*sizeof(double));
    cblas_dgemv(CblasRowMajor,CblasNoTrans,m*N,n,1,controller->K,n,x0,1,1,uopt,1);
  }
  print_cpuinfo("cblas_dgemv(K,x0)", t, NTESTS);
  struct timespec clk[NTESTS];
  for (uint32_t i=0; i<NTESTS; i++){
    clock_gettime(CLOCK_MONOTONIC, &clk[i]);
    memset(uopt,0x00,GPQHE_HORIZON*GPQHE_DIMU*sizeof(double));
    cblas_dgemv(CblasRowMajor,CblasNoTrans,m*N,n,1,controller->K,n,x0,1,1,uopt,1);
  }
  print_clkinfo("cblas_dgemv(K,x0)", clk, NTESTS);
  uint64_t maxrss[NTESTS];
  for (uint32_t i=0; i<NTESTS; i++){
    maxrss[i] = rdmaxrss();
    memset(uopt,0x00,GPQHE_HORIZON*GPQHE_DIMU*sizeof(double));
    cblas_dgemv(CblasRowMajor,CblasNoTrans,m*N,n,1,controller->K,n,x0,1,1,uopt,1);
  }
  print_meminfo("cblas_dgemv(K,x0)", maxrss, NTESTS);
#endif
  /* state evolution */
  sprintf(name,
    "encrypted_mpc/finite_horizon_optimal_input_[%g,%g]_plain",
    x0[0],x0[1]);
  actuate(name, x0, uopt, controller);
  /* extend K */
  size_t d_half=GPQHE_N/2;
  COMPLEX K_ext[(GPQHE_N/2)*(GPQHE_N/2)] = {0};
  for (size_t i=0; i<m*N*n; i++)
    K_ext[(i/n)*d_half+i%n]=controller->K[i];
  /* encode & encrypt */
  COMPLEX msg[GPQHE_N/2] = {0};
  for (size_t i=0; i<GPQHE_DIMX; i++)
    msg[i]=x0[i];
  ckks_pt pt;
  ckks_ct cx,cu;
  he_ckks_encode(&pt, d_half, msg, ctx);
  he_ckks_encrypt_pk(&cx, &pt, pk, ctx);
  /* cuopt */
  double cuopt[GPQHE_HORIZON*GPQHE_DIMU] = {0}; /* N*m */
#if GPQHE_PMU
  for (uint32_t i=0; i<NTESTS; i++){
    t[i] = rdtsc();
    he_ckks_gemv(&cu, K_ext, &cx, evk->rot_key, ctx);
  }
  print_cpuinfo("he_ckks_gemv(K,Enc(x0))", t, NTESTS);
  for (uint32_t i=0; i<NTESTS; i++){
    clock_gettime(CLOCK_MONOTONIC, &clk[i]);
    he_ckks_gemv(&cu, K_ext, &cx, evk->rot_key, ctx);
  }
  print_clkinfo("he_ckks_gemv(K,Enc(x0))", clk, NTESTS);
  for (uint32_t i=0; i<NTESTS; i++){
    maxrss[i] = rdmaxrss();
    he_ckks_gemv(&cu, K_ext, &cx, evk->rot_key, ctx);
  }
  print_meminfo("he_ckks_gemv(K,Enc(x0))", maxrss, NTESTS);
#endif
  he_ckks_gemv(&cu, K_ext, &cx, evk->rot_key, ctx);
  he_ckks_decrypt(&pt, &cu, sk, ctx);
  he_ckks_decode(d_half, msg, &pt, ctx);
  for (uint32_t i=0; i<m*N; i++)
  {
    assert(cimag(msg[i])<FLT_EPSILON);
    cuopt[i]=creal(msg[i]);
  }
  /* state evolution */
  sprintf(name,
    "encrypted_mpc/finite_horizon_optimal_input_[%g,%g]_enc_sk",
    x0[0],x0[1]);
  actuate(name, x0, cuopt, controller);
  /* deviation */
  printf("sk encrypted: log10(||gemv(K,x0)-gemv(K,Enc(x0))||) = %.15g\n",
    log10(blas_dnrmmax_dist(GPQHE_HORIZON*GPQHE_DIMU, uopt, cuopt)));
  TEST_END();
}
#endif
int main(int argc __attribute__((unused)),
         const char *argv[] __attribute__((unused)))
{
#if 0
  struct controller_ctx controller;
#if GPQHE_PRECOMP
  linear_controller_ctx_build(&controller);
  linear_controller_ctx_write(&controller);
#endif
  linear_controller_ctx_load(&controller);
#endif

  gcry_control (GCRYCTL_DISABLE_SECMEM, 0);
  gcry_control (GCRYCTL_ENABLE_QUICK_RANDOM, 0);
  gcry_control (GCRYCTL_INITIALIZATION_FINISHED, 0);
  #if 0
  test_convert();
  test_bit_operations();
  test_blas();
  test_crt();
  test_ntt();
  test_fft();
  test_embedding();
  test_poly_mod();
  test_poly_add();
  test_poly_sub();
  test_poly_scal();
  test_poly_rot();
  //test_polynomial();
  test_poly_multiply();
  #endif

  /* below requires initialize context */
  
  he_ckks_ctx_build(0);

#if 1
  test_encode();
  he_ckks_ctx_show_params();
#endif

  int8_t sk[GPQHE_INDCPA_SECRETKEYBYTES]={0};
  struct ckks_pk pk;
  struct ckks_pk rlk __attribute__((unused));
  struct ckks_pk ck  __attribute__((unused));
  struct ckks_pk rk[GPQHE_CKKS_SLOT] __attribute__((unused));
  indcpa_keypair(sk, &pk);

#if 1
  test_encrypt_decrypt(sk, &pk);
  test_addsub(sk, &pk);
#endif

#if 1
  he_ckks_genrlk(&rlk, sk);
  test_multiply(sk, &pk, &rlk);
  test_exp(sk, &pk, &rlk);
#endif

#if 1
  he_ckks_genck(&ck, sk);
  test_conj(sk, &pk, &ck);
#endif

#if 0
  he_ckks_genrk(rk, sk);
  test_rot(sk, &pk, rk);
#endif

#if 1
  he_ckks_genrk(rk, sk);
  test_gemv(sk, &pk, rk);
#endif

#if 0

#ifndef GPQHE_GENRLK
  he_ckks_genrlk(&rlk, sk);
  #define GPQHE_GENRLK
#endif

#ifndef GPQHE_GENCK
  he_ckks_genck(&ck, sk);
  #define GPQHE_GENCK
#endif

  test_rlsin(sk, &pk, &rlk, &ck);
#endif

#if 0

#ifndef GPQHE_GENCK
  he_ckks_genck(&ck, sk);
  #define GPQHE_GENCK
#endif

#ifndef GPQHE_GENRK
  he_ckks_genrk(rk, sk);
  #define GPQHE_GENRK
#endif

  test_coeff_to_slot(sk, &pk, &ck, rk);
  test_slot_to_coeff(sk, &pk, &ck, rk);
#endif

#if 0

#ifndef GPQHE_GENRLK
  he_ckks_genrlk(&rlk, sk);
  #define GPQHE_GENRLK
#endif

#ifndef GPQHE_GENCK
  he_ckks_genck(&ck, sk);
  #define GPQHE_GENCK
#endif

#ifndef GPQHE_GENRK
  he_ckks_genrk(rk, sk);
  #define GPQHE_GENRK
#endif

  test_bootstrap(sk, &pk, &rlk, &ck, rk);
#endif

#if 0

#ifndef GPQHE_GENRLK
  he_ckks_genrlk(&rlk, sk);
  #define GPQHE_GENRLK
#endif

  test_inv(sk, &pk, &rlk);
#endif


#if 1

#ifndef GPQHE_GENRLK
  he_ckks_genrlk(&rlk, sk);
  #define GPQHE_GENRLK
#endif

  test_sqrt(sk, &pk, &rlk);
#endif

#if 0
  test_he_ckks_sin(&sk, &pk, &evk, &ctx);
  test_finite_horizon_optimal_input(argc, argv, &sk, &pk, &evk, &ctx, &controller);
#endif
  return 0;
}
