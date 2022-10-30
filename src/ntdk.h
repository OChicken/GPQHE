/*
 * Number Theory Development Kit based on Libgcrypt
 * Copyright (C) shouran.ma@rwth-aachen.de
 *
 * This file is part of Kyfher.
 *
 * Kyfher is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * Kyfher is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this program; if not, see <http://www.gnu.org/licenses/>.
 */

/* begin ntdk.h */

#ifndef __NTDK_H__
#define __NTDK_H__

#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>

/* math */
#ifndef __USE_XOPEN
  #define __USE_XOPEN
#endif

#include <complex.h>
typedef double _Complex COMPLEX;
typedef struct gcry_mpi *MPI;

#include <math.h>
#include "config.h"
#define GPQHE_N 128//16384
#define GPQHE_L 20
#define GPQHE_P (1ULL<<50)
#define GPQHE_HW 64
#define GPQHE_SIGMA (double)3.2
#define GPQHE_CKKS_ITER 10

/* libgcrypt */
#include <gcrypt.h>

#ifndef KYBER_K
#define KYBER_K 3 /* Change this for different security strengths */
#endif

/* polynomial degree
 * The following lines are defined in Kyber/Reference_Implementation/crypto_kem/kyber1024:
 * params.h:37:#define KYBER_SYMBYTES 32   // size in bytes of hashes, and seeds
 * params.h:59:#define KYBER_INDCPA_MSGBYTES       KYBER_SYMBYTES
 * poly.c:147:void poly_frommsg(poly *r, const uint8_t msg[KYBER_INDCPA_MSGBYTES])
 * poly.c:152:#if (KYBER_INDCPA_MSGBYTES != KYBER_N/8)
 * poly.c:153:#error "KYBER_INDCPA_MSGBYTES must be equal to KYBER_N/8 bytes!"
 * poly.c:154:#endif
 * In practice, the polynomial degree is set to be 1024.
 */
#define KYBER_N 256

#if KYBER_K == 2
#define KYBER_ETA1 3
#define KYBER_POLYCOMPRESSEDBYTES    128
#define KYBER_POLYVECCOMPRESSEDBYTES (KYBER_K * 320)
#elif KYBER_K == 3
#define KYBER_ETA1 2
#define KYBER_POLYCOMPRESSEDBYTES    128
#define KYBER_POLYVECCOMPRESSEDBYTES (KYBER_K * 320)
#elif KYBER_K == 4
#define KYBER_ETA1 2
#define KYBER_POLYCOMPRESSEDBYTES    160
#define KYBER_POLYVECCOMPRESSEDBYTES (KYBER_K * 352)
#endif

#define NUM_TAYLOR_ITER (size_t)7    /* number of Taylor iteration */
#define PBND 59
#define LOGQ 800
#define USE_RNS 1
#define KNOWN_GENERATOR 1

BEGIN_DECLS

/*******************************************************************************
 * Type Conversion Utils                                                       *
 ******************************************************************************/

/* Show sexp */
void show_mpi (MPI a);

/* sprint hex */
char *sprint_hex(MPI N);

/* Convert MPI to int64_t */
int64_t gcry_mpi_to_int64 (MPI a);

/* Convert MPI to uint64_t */
uint64_t gcry_mpi_to_uint64 (MPI a);

/* Convert MPI to int128_t */
__int128_t gcry_mpi_to_int128 (MPI a);

/* Convert MPI to uint128_t */
__uint128_t gcry_mpi_to_uint128 (MPI a);

/* Convert int64_t to MPI */
MPI gcry_int64_to_mpi (int64_t a);

/* Convert int128_t to MPI */
MPI gcry_int128_to_mpi (__int128_t a);

/* buffer to unsigned long */
uint64_t mpi_buf2uint64 (const char *buf, const size_t len);

/* buffer to unsigned long */
uint64_t ntdk_buf_to_u64 (const unsigned char *buf, const size_t len);

/*******************************************************************************
 * Arithmetic Utils                                                            *
 ******************************************************************************/

/* MPI square - bisection method */
MPI gcry_mpi_sqrt(MPI x);

/* Return floor(log2(x)). */
size_t ntdk_log2(uint64_t x);

/* Get the number of bits required to represent an unsigned 64-bit number x */
size_t ntdk_get_nbits(uint64_t x);

/* Bit reversal operation, with width can be specified. */
long mpi_br_width (long x, size_t w);

/* Bit reversal operation to an unsigned 32-bit number */
uint32_t ntdk_br_u32(uint32_t x);

/* Bit reversal operation to an unsigned 64-bit number */
uint64_t ntdk_br_u64(uint64_t x);

/* square root to unsigned 64-bit number */
uint64_t ntdk_sqrt_u64(uint64_t x);

/* power */
uint64_t ntdk_pow(uint64_t b, uint64_t e);

/* inverse in unsigned 64-bit */
uint64_t ntdk_inv_u64(uint64_t x);

/* multiply mod to unsigned 64-bit inputs */
uint64_t ntdk_mulm(uint64_t a, uint64_t b, uint64_t m);

/* multiply mod to unsigned 64-bit inputs */
uint64_t ntdk_mulm_barrett(uint64_t a, uint64_t b, uint64_t p, uint64_t pr);

/* power mod */
uint64_t ntdk_powm(uint64_t b, uint64_t e, uint64_t m);

/* signed mod */
void ntdk_smod(MPI *a, const MPI m);

/* rounding to the nearest integer division */
void ntdk_div_round(MPI q, MPI dividend, MPI divisor);

/* M-th root of unity */
uint64_t ntdk_m_root_of_unity(uint64_t m, uint64_t p);

/*******************************************************************************
 * Polynomial                                                                  *
 ******************************************************************************/

/* Polynomials
 * Elements of R_q = Z_q[X]/(X^n + 1). Represents polynomial
 * coeffs[0] + X*coeffs[1] + X^2*xoeffs[2] + ... + X^{n-1}*coeffs[n-1] */
/* polynomial of integer coefficients (no suffix by default) */
typedef struct{
  size_t d;
  size_t m; /* cyclotomic index */
  uint64_t *p;
  uint64_t *pr;
  uint64_t *p_inv;
  uint64_t *d_inv_scaled;
  uint64_t **rootpow_scaled;
  uint64_t **rootpow_inv_scaled;
  MPI *p_prod;
  MPI *ph_prod;
  MPI *ph;
  MPI *ph_inv_modp;
} poly_ctx;

void poly_ctx_build(poly_ctx *poly, size_t d);

//void polyz_build(ntt_ctx *ntt, size_t d, MPI p, MPI g);

/* polynomial of real coefficients (suffix "d" implies "double") */
typedef struct{
  size_t d;
  double *coeffs;
} polyr;


/*******************************************************************************
 * BLAS                                                                        *
 * ****************************************************************************/

/* matrix vector multiply */
void blas_lgemv(size_t m, size_t n, long *y, long *A, long *x);
void blas_zgemv(size_t m, size_t n, _Complex double *y, _Complex double *A, _Complex double *x);

/* vector addition */
void blas_xpy(size_t m, long *r, long *x, long *y);

/* complex double vector addition */
void blas_zxpy(size_t m, _Complex double *r, _Complex double *x, _Complex double *y);

/* complex double vector subtraction */
void blas_zxsy(size_t m, _Complex double *r, _Complex double *x, _Complex double *y);

/* vector multiplication */
void blas_xmy(size_t m, long *r, long *x, long *y);

/* complex double vector multiplication */
void blas_zxmy(size_t m, _Complex double *r, _Complex double *x, _Complex double *y);

/* scalar multiply */
void blas_scal(size_t m, long *y, long alpha, long *x);

/* fetch diagonal of the matrix A
 * A: 2D array of length m*m
 * diag: 1D array of length m     */
void blas_ldiag(size_t m, long *diag, size_t index, long *A);
void blas_ddiag(size_t m, double *diag, size_t index, double *A);
void blas_zdiag(size_t m, _Complex double *diag, size_t index, _Complex double *A);

/* rotate
 * x and y: 1D array 
 * m: length of x and y */
void blas_lrot(size_t m, long *y, int64_t rot, long *x);
void blas_drot(size_t m, double *y, int64_t rot, double *x);
void blas_zrot(size_t m, _Complex double *y, int64_t rot, _Complex double *x);

/* conjugate */
void blas_conj(size_t m, _Complex double *y, _Complex double *x);

/* transpose */
void blas_trans(size_t m, size_t n, long *AT, long *A);

/* Return 1 if two MPI type vectors are equal; 0 if else. */
uint8_t blas_mpivec_equal(size_t m, MPI *y, MPI *x);

/* Return 1 if two long type vectors are equal; 0 if else. */
uint8_t blas_lvec_equal(size_t m, long *y, long *x);

/* Return 1 if two complex double type vectors are equal (below the threshold); 0 if else. */
uint8_t blas_zvec_equal(size_t m, _Complex double *y, _Complex double *x, double threshold);

/* Infinity norm of a _Complex double type vector */
double blas_dznrmmax(size_t m, _Complex double *z);

/* Infinity norm of the difference of two _Complex double type vectors */
double blas_dznrmmax_dist(size_t m, _Complex double *y, _Complex double *x);

/* Infinity norm of the difference of two double type vectors */
double blas_dnrmmax_dist(size_t m, double *y, double *x);

/*******************************************************************************
 * Sample and Statistics                                                       *
 ******************************************************************************/

/* sample uniform complex random vector with each coefficient in [0,1] */
void sample_z01vec(size_t m, _Complex double *vec);

/* sample centered binomial distribution 
 * Given an array of uniformly random bytes, compute polynomial with 
 * coefficients distributed according to a centered binomial distribution with
 * parameter eta=2 */
#if KYBER_ETA1 == 2
void sample_cbd2_eta1(size_t m, int16_t *vec);
#endif

/* sample discrete gaussian */
void sample_discrete_gaussian(size_t m, int16_t *vec, double sigma);

/* sample hamming weight vector, whose hamming weight is exactly h */
void sample_hwt(size_t m, int8_t *vec, size_t h);

/* sample zero-centered vector, */
void sample_zero_center(size_t m, int16_t *vec);

/* sample uniform_u64*/
void sample_uniform_u64(size_t m, uint64_t *r, uint64_t q);

/*******************************************************************************
 * CRT                                                                         *
 ******************************************************************************/

/* CRT context */
typedef struct{
  size_t L;       /* length of ps, ms, ms_inv */
  MPI *p;  /* prime list */
  MPI *g;  /* generator list */
  MPI *m;  /* list of mi, where mi = M/pi */
  MPI *c;  /* ci = mi*inv(mi, M) */
  MPI M;   /* products of all primes in prime list */
} crt_ctx;

/* Build CRT context */
void crt_build(crt_ctx *crt, size_t size_p, size_t num_p);

/* transform a number `n` to CRT representation list `a`
 * Here, `n` is the coefficient of a polynomial, so signed */
void crt_fwd(MPI *a, MPI n, crt_ctx *crt);

/* reconstruct `n` from vectors `a` in CRT rep */
void crt_inv(MPI *n, MPI *a, crt_ctx *crt);

/* NTT (Number Theoretic Transformation) context */
typedef struct{
  size_t d;
  MPI p;
  MPI g;
  MPI *zetas;
  MPI *zetas_inv;
} ntt_ctx;

/* Build NTT context with PRIME modulus p
 * g is the generator of p. If the exact value of g is not know, set g=0. */
void ntt_build(ntt_ctx *ntt, size_t d, MPI p, MPI g);

/* forward NTT/FTT (Fermat Theoretic Transform) (modulus p must be prime) */
void ntt_fwd(size_t m, MPI *X, const MPI *x, const MPI p, const ntt_ctx *ntt);

/* inverse NTT/FTT (Fermat Theoretic Transform) (modulus must be prime) */
void ntt_inv(size_t m, MPI *x, const MPI *X, const MPI p, const ntt_ctx *ntt);

/* FFT context */
typedef struct{
  size_t L;
  _Complex double *zetas;
  int64_t *rot_group;
} fft_ctx;

/* Build FFT context */
void fft_build(fft_ctx *fft, size_t d);

/* forward FFT */
void fft_fwd(size_t m, _Complex double *X, const _Complex double *x, const fft_ctx *fft);

/* inverse FFT */
void fft_inv(size_t m, _Complex double *x, const _Complex double *X, const fft_ctx *fft);

/* forward canonical embedding */
void emb_fwd(size_t m, _Complex double *X, const _Complex double *x, const fft_ctx *fft);

/* inverse canonical embedding */
void emb_inv(size_t m, _Complex double *x, const _Complex double *X, const fft_ctx *fft);

/* Polynomials
 * Elements of R_q = Z_q[X]/(X^n + 1). Represents polynomial
 * coeffs[0] + X*coeffs[1] + X^2*xoeffs[2] + ... + X^{n-1}*coeffs[n-1] */
/* polynomial of integer coefficients (no suffix by default) */
typedef struct{
  size_t d;
  int64_t *coeffs;
} poly;

/* polynomial of real coefficients (suffix "d" implies "double") */
typedef struct{
  size_t d;
  double *coeffs;
} polyd;

/* poly modulus */
void poly_mod(MPI *r, const MPI Q);

/* poly signed modulus: mod to the range [-q/2, q/2); a.k.a. small mod */
void poly_smod(MPI *r, const MPI Q);

/* poly round off to the nearest integer
 * Rounds all the x's coefficients to the nearest integer, where |x| = n + 0.5
 * rounds to |x| = n (i.e. 0.5 rounds to 0 and -1.5 rounds to -1). */
void poly_round(size_t d, MPI *r, const polyd *x);

/* poly addition
 * if Q!=NULL and flag==0, r[i] =  mod(x[i]+y[i], Q);
 * if Q!=NULL and flag==1, r[i] = smod(x[i]+y[i], Q);
 * if Q==NULL              r[i] =      x[i]+y[i]      */
void poly_xpy(MPI *r,
  const MPI *x, const MPI *y,
  const MPI Q, unsigned int flag);

/* poly subtraction
 * if Q!=NULL and flag==0, r[i] =  mod(x[i]-y[i], Q);
 * if Q!=NULL and flag==1, r[i] = smod(x[i]-y[i], Q);
 * if Q==NULL              r[i] =      x[i]-y[i]      */
void poly_xsy(MPI *r,
  const MPI *x, const MPI *y,
  const MPI Q, unsigned int flag);

/* poly multiplication
 * Multiplies the poly a and b inside the ring R_a naively in O(n^2) time.
 * if Q!=NULL and flag==0, r[i] =  mod(x[i]*y[i], Q);
 * if Q!=NULL and flag==1, r[i] = smod(x[i]*y[i], Q);
 * if Q==NULL              r[i] =      x[i]*y[i]      */
void poly_xmy(MPI *r,
  const MPI *x, const MPI *y,
  const MPI Q, unsigned int flag);

/* poly multiplication with CRT
 * Multiplies the poly a and b inside the ring by splitting it into CRT subrings
 * for the primes given. For each subring, we multiply using NTT and recombine
 * with CRT. The multiplications costs O(n^2) time. */
void poly_xmy_crt(size_t d, MPI *r,
  const MPI *x, const MPI *y, crt_ctx *crt, ntt_ctx *ntts);

/* poly multiplication with FFT */
void poly_xmy_fft(polyd *r, const polyd *x, const polyd *y);

/* poly scalar multiply */
void poly_scal(size_t d, MPI *r, const long alpha);

/* poly rotation */
void poly_rot(MPI *r, const MPI *x, size_t rot);

/* poly conjugation */
void poly_conj(MPI *r, const MPI *x);

/* poly evaluation */
MPI poly_eval(size_t d, MPI *r, const int64_t x0);

/* Determine Euler totient and the factorization of the input number N.

   Args:
     gcrypt_mpi_t N: the number to be factorized.
     gcrypt_mpi_t *phi: Euler totient phi(N).
     gcrypt_mpi_t **factors: store the list of factorize(N).
 
   Return:
     gcrypt_error_t
 */
gcry_error_t
gcry_ntdk_totient(MPI N,
                  MPI *phi,
                  MPI **factors);

/* Determine the group generator / primitive root of the input number N.

   Args:
     gcrypt_mpi_t *r_g: the generator.
     gcrypt_mpi_t N: the modulus.
 
   Return:
     gcrypt_error_t
 */
gcry_error_t
gcry_ntdk_group_generator (MPI *r_g, MPI N);

END_DECLS

#endif /* __NTDK_H__ */

/* end ntdk.h */
