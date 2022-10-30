/*
 * kem - key encapsulation mechanism.
 * Copyright (C) shouran.ma@rwth-aachen.de
 *
 * This file is part of GPQHE.
 *
 * GPQHE is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * GPQHE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this program; if not, see <http://www.gnu.org/licenses/>.
 */

#include "fhe.h"

BEGIN_DECLS

/* poly.h */
extern struct poly_ctx polyctx;

/* fhe.h */
extern struct he_ctx hectx;

/* types.c */
extern void s64_to_mpi(MPI r, int64_t a);
extern void mpi_smod(MPI r, const MPI q, const MPI qh);

/* rng.c */
extern void randombytes(uint8_t *x,size_t xlen);

/* sample.h */
extern void sample_hwt(int64_t *vec, const size_t m);
extern void sample_discrete_gaussian(int64_t *vec, const size_t m);
extern void sample_zero_center(int64_t *vec, const size_t m);
extern void sample_uniform(poly_mpi_t *r, const MPI q);

/** he_keypair - sk and pk key pair */
void he_keypair(he_pk_t *pk, poly_mpi_t *sk)
{
  /* local variables */
  MPI qL = mpi_copy(hectx.q[hectx.L]);
  MPI qh = mpi_copy(hectx.qh[hectx.L]);
  unsigned int n = polyctx.n;
  unsigned int dim = mpi_get_nbits(qL)/GPQHE_LOGP+1;
  printf("Generating sk and pk ... ");
  fflush(stdout);
  /* secret key */
  int64_t sbuf[n];
  memset(sbuf, 0, sizeof(sbuf));
  sample_hwt(sbuf, n);
  for (size_t i=0; i<n; i++)
    s64_to_mpi(sk->coeffs[i], sbuf[i]);
  /* sample error for public key */
  MPI e = mpi_new(0);
  int64_t ebuf[n];
  memset(ebuf, 0, sizeof(ebuf));
  sample_discrete_gaussian(ebuf, n);
  /* main */
  sample_uniform(&pk->p1, qL);
  poly_mul(&pk->p0, sk, &pk->p1, dim, qL);
  for (unsigned int i=0; i<polyctx.n; i++) {
    s64_to_mpi(e, ebuf[i]);
    mpi_neg(pk->p0.coeffs[i], pk->p0.coeffs[i]);
    mpi_addm(pk->p0.coeffs[i], pk->p0.coeffs[i], e, qL);
    mpi_smod(pk->p0.coeffs[i], qL, qh);
    mpi_smod(pk->p1.coeffs[i], qL, qh); /* add this or not does not affect the correctness */
  }
  /* release */
  mpi_release(e);
  mpi_release(qL);
  mpi_release(qh);
  printf("done.\n");
}

/** switching keys */
static void he_genswk(he_pk_t *swk, const poly_mpi_t *sp, const poly_mpi_t *sk)
{
  /* local variables */
  MPI P    = mpi_copy(hectx.P);
  MPI PqL  = mpi_new(0);
  MPI PqLh = mpi_new(0);
  mpi_mul(PqL,  P, hectx.q [hectx.L]);
  mpi_mul(PqLh, P, hectx.qh[hectx.L]);
  unsigned int n = polyctx.n;
  unsigned int dim = mpi_get_nbits(PqL)/GPQHE_LOGP+1;
  /* error */
  MPI e = mpi_new(0);
  int64_t ebuf[n];
  sample_discrete_gaussian(ebuf, n);
  /* main */
  for (size_t i=0; i<n; i++)
    mpi_mul(sp->coeffs[i], sp->coeffs[i], P);
  sample_uniform(&swk->p1, PqL);
  poly_mul(&swk->p0, &swk->p1, sk, dim, PqL);
  for (unsigned int i=0; i<n; i++) {
    s64_to_mpi(e, ebuf[i]);
    mpi_neg(swk->p0.coeffs[i], swk->p0.coeffs[i]);                      /* swk.p0 = -swk.p1*sk */
    mpi_add(swk->p0.coeffs[i], swk->p0.coeffs[i], e);                   /* swk.p0 = -swk.p1*sk+e */
    mpi_addm(swk->p0.coeffs[i], swk->p0.coeffs[i], sp->coeffs[i], PqL); /* swk.p0 = -swk.p1*sk+e+Psp mod PqL */
    mpi_smod(swk->p0.coeffs[i], PqL, PqLh);
    mpi_smod(swk->p1.coeffs[i], PqL, PqLh);
  }
  mpi_release(P);
  mpi_release(PqL);
  mpi_release(PqLh);
  mpi_release(e);
}

void he_genrlk(he_pk_t *rlk, const poly_mpi_t *sk)
{
  /* local variables */
  MPI q = mpi_copy(hectx.q[hectx.L]);
  poly_mpi_t s2;
  poly_mpi_alloc(&s2);
  /* main */
  printf("Generating rlk ... ");
  fflush(stdout);
  /* relinearization key */
  poly_mul(&s2, sk, sk, mpi_get_nbits(q)/GPQHE_LOGP+1, q);
  he_genswk(rlk, &s2, sk);
  printf("done.\n");
  /* release */
  mpi_release(q);
  poly_mpi_free(&s2);
}

/** conjugate key */
void he_genck(he_pk_t *ck, const poly_mpi_t *sk)
{
  /* local variables */
  poly_mpi_t ck_;
  poly_mpi_alloc(&ck_);
  /* conjugate key */
  printf("Generating ck ... ");
  fflush(stdout);
  poly_conj(&ck_, sk);
  he_genswk(ck, &ck_, sk);
  printf("done.\n");
  /* release */
  poly_mpi_free(&ck_);
}

void he_genrk(he_pk_t *rk, const poly_mpi_t *sk)
{
  /* local variables */
  poly_mpi_t rk_;
  poly_mpi_alloc(&rk_);
  /* rotation keys */
  printf("Generating rk ... ");
  fflush(stdout);
  for (unsigned int rot=0; rot<hectx.slots; rot++) {
    poly_rot(&rk_, sk, rot);
    he_genswk(&rk[rot], &rk_, sk);
  }
  printf("done.\n");
  /* release */
  poly_mpi_free(&rk_);
}

END_DECLS
