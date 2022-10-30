/*
 * FHE operation collections of multiplication and rescale.
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
#include <math.h>

BEGIN_DECLS

/* poly.h */
extern struct poly_ctx polyctx;

/* fhe.h */
extern struct he_ctx hectx;

/* types.c */
extern void s64_to_mpi(MPI r, int64_t a);
extern void mpi_smod(MPI r, const MPI q, const MPI qh);

/* ntt.c */
extern void    ntt(uint64_t a[], const struct rns_ctx *rns);
extern void invntt(uint64_t a[], const struct rns_ctx *rns);

/* rns.c */
extern void rns_decompose(uint64_t ahat[], const MPI a[], const struct rns_ctx *rns);
extern void rns_reconstruct(MPI a[], const uint64_t ahat[], const unsigned int i, const struct rns_ctx *rns);

/** Relinearizes a dim3 ciphertext back down to dim2 */
static void he_relin(struct he_ct *ct,
  const poly_mpi_t *d0, const poly_mpi_t *d1, const poly_mpi_t *d2,
  const he_pk_t *rlk, const MPI ql)
{
  /* local variables */
  MPI qlh = mpi_copy(hectx.qh[ct->l]);
  MPI P   = mpi_copy(hectx.P);
  MPI Pql = mpi_new(0);
  mpi_mul(Pql, P, ql);
  unsigned int n = polyctx.n;
  unsigned int dim = (mpi_get_nbits(ql)+mpi_get_nbits(Pql))/GPQHE_LOGP+1;
  /* main */
  poly_mul(&ct->c0, d2, &rlk->p0, dim, Pql);
  poly_mul(&ct->c1, d2, &rlk->p1, dim, Pql);
  for (size_t i=0; i<n; i++) {
    mpi_rdiv(ct->c0.coeffs[i], ct->c0.coeffs[i], P);
    mpi_rdiv(ct->c1.coeffs[i], ct->c1.coeffs[i], P);
    mpi_mod(ct->c0.coeffs[i], ct->c0.coeffs[i], ql);
    mpi_mod(ct->c1.coeffs[i], ct->c1.coeffs[i], ql);
    mpi_addm(ct->c0.coeffs[i], ct->c0.coeffs[i], d0->coeffs[i], ql);
    mpi_addm(ct->c1.coeffs[i], ct->c1.coeffs[i], d1->coeffs[i], ql);
    mpi_smod(ct->c0.coeffs[i], ql, qlh);
    mpi_smod(ct->c1.coeffs[i], ql, qlh);
  }
  /* release */
  mpi_release(qlh);
  mpi_release(P);
  mpi_release(Pql);
}

/** ct_dest = ct1*ct2 */
void he_mul(he_ct_t *ct, const he_ct_t *ct1, const he_ct_t *ct2, const he_pk_t *rlk)
{
  assert(ct1->l==ct2->l);
  /* ct init */
  ct->l  = ct1->l;
  ct->nu = ct1->nu*ct2->nu;
  ct->B  = (ct1->nu)*(ct2->B) + (ct2->nu)*(ct1->B) + (ct1->B)*(ct2->B)
         + hectx.bnd.Bmult[ct->l];
  /* local variables */
  MPI q  = mpi_copy(hectx.q [ct->l]);
  unsigned int n = polyctx.n;
  unsigned int dim = (mpi_get_nbits(q)*2)/GPQHE_LOGP+1;
  /* alloc */
  poly_mpi_t d0, d1, d2;
  poly_mpi_alloc(&d0);
  poly_mpi_alloc(&d1);
  poly_mpi_alloc(&d2);
  poly_rns_t ct1c0hat, ct1c1hat, ct2c0hat, ct2c1hat;
  poly_rns_alloc(&ct1c0hat, 1);
  poly_rns_alloc(&ct1c1hat, 1);
  poly_rns_alloc(&ct2c0hat, 1);
  poly_rns_alloc(&ct2c1hat, 1);
  poly_rns_t d0hat, d1hat, d2hat, d11hat, d12hat;
  poly_rns_alloc(&d0hat, dim);
  poly_rns_alloc(&d1hat, dim);
  poly_rns_alloc(&d2hat, dim);
  poly_rns_alloc(&d11hat, 1);
  poly_rns_alloc(&d12hat, 1);
  /* main - cross terms */
  struct rns_ctx *rns = polyctx.rns;
  for (unsigned int d=0; d<dim; d++) {
    rns_decompose(ct1c0hat.coeffs, ct1->c0.coeffs, rns);
    rns_decompose(ct1c1hat.coeffs, ct1->c1.coeffs, rns);
    rns_decompose(ct2c0hat.coeffs, ct2->c0.coeffs, rns);
    rns_decompose(ct2c1hat.coeffs, ct2->c1.coeffs, rns);
    ntt(ct1c0hat.coeffs, rns);
    ntt(ct1c1hat.coeffs, rns);
    ntt(ct2c0hat.coeffs, rns);
    ntt(ct2c1hat.coeffs, rns);
    poly_rns_mul( &d0hat.coeffs[d*n], ct1c0hat.coeffs, ct2c0hat.coeffs, rns);
    invntt(&d0hat.coeffs[d*n], rns); /* d0 = ct1.c0*ct2.c0 */
    poly_rns_mul( &d2hat.coeffs[d*n], ct1c1hat.coeffs, ct2c1hat.coeffs, rns);
    invntt(&d2hat.coeffs[d*n], rns); /* d2 = ct1.c1*ct2.c1 */
    poly_rns_mul(d11hat.coeffs, ct1c0hat.coeffs, ct2c1hat.coeffs, rns);
    invntt(d11hat.coeffs, rns); /* d11 = ct1.c0*ct2.c1 */
    poly_rns_mul(d12hat.coeffs, ct1c1hat.coeffs, ct2c0hat.coeffs, rns);
    invntt(d12hat.coeffs, rns); /* d12 = ct1.c1*ct2.c0 */
    poly_rns_add(&d1hat.coeffs[d*n], d11hat.coeffs, d12hat.coeffs, rns);
    if (d<dim-1)
      rns=rns->next;
  }
  poly_rns2mpi(&d0, &d0hat, rns, q); /* d0 = ct1.c0*ct2.c0 */
  poly_rns2mpi(&d1, &d1hat, rns, q); /* d1 = ct1.c0*ct2.c1 + ct1.c1*ct2.c0 */
  poly_rns2mpi(&d2, &d2hat, rns, q); /* d2 = ct1.c1*ct2.c1 */
  /* main - relinearize */
  he_relin(ct, &d0, &d1, &d2, rlk, q);
  /* release */
  mpi_release(q);
  poly_rns_free(&ct1c0hat);
  poly_rns_free(&ct1c1hat);
  poly_rns_free(&ct2c0hat);
  poly_rns_free(&ct2c1hat);
  poly_rns_free(&d0hat);
  poly_rns_free(&d1hat);
  poly_rns_free(&d2hat);
  poly_rns_free(&d11hat);
  poly_rns_free(&d12hat);
  poly_mpi_free(&d0);
  poly_mpi_free(&d1);
  poly_mpi_free(&d2);
}

/** ct_dest = ct*pt */
void he_mulpt(he_ct_t *dest, const he_ct_t *src, const he_pt_t *pt)
{
  /* ct init */
  dest->l  = src->l;
  dest->nu = src->nu*pt->nu;
  dest->B  = src->B*pt->nu;
  /* local varables */
  MPI q  = mpi_copy(hectx.q [dest->l]);
  unsigned int n = polyctx.n;
  unsigned int dim = (mpi_get_nbits(q)+log2(pt->nu))/GPQHE_LOGP+1;
  /* alloc */
  poly_rns_t pthat, srcc0hat, srcc1hat, destc0hat, destc1hat;
  poly_rns_alloc(&pthat, 1);
  poly_rns_alloc(&srcc0hat, 1);
  poly_rns_alloc(&srcc1hat, 1);
  poly_rns_alloc(&destc0hat, dim);
  poly_rns_alloc(&destc1hat, dim);
  /* main */
  struct rns_ctx *rns = polyctx.rns;
  for (unsigned int d=0; d<dim; d++) {
    rns_decompose(pthat.coeffs, pt->m.coeffs, rns);
    rns_decompose(srcc0hat.coeffs, src->c0.coeffs, rns);
    rns_decompose(srcc1hat.coeffs, src->c1.coeffs, rns);
    ntt(pthat.coeffs, rns);
    ntt(srcc0hat.coeffs, rns);
    ntt(srcc1hat.coeffs, rns);
    poly_rns_mul(&destc0hat.coeffs[d*n], srcc0hat.coeffs, pthat.coeffs, rns);
    invntt(&destc0hat.coeffs[d*n], rns);
    poly_rns_mul(&destc1hat.coeffs[d*n], srcc1hat.coeffs, pthat.coeffs, rns);
    invntt(&destc1hat.coeffs[d*n], rns);
    if (d<dim-1)
      rns=rns->next;
  }
  poly_rns2mpi(&dest->c0, &destc0hat, rns, q);
  poly_rns2mpi(&dest->c1, &destc1hat, rns, q);
  /* release */
  mpi_release(q);
  poly_rns_free(&pthat);
  poly_rns_free(&srcc0hat);
  poly_rns_free(&srcc1hat);
  poly_rns_free(&destc0hat);
  poly_rns_free(&destc1hat);
}

/** rescale */
void he_rs(struct he_ct *ct)
{
  /* ct init */
  ct->l -= 1;
  ct->nu /= hectx.Delta;
  ct->B = ct->B/hectx.Delta + hectx.bnd.Brs;
  /* local varables */
  unsigned int n = polyctx.n;
  MPI q  = mpi_copy(hectx.q [ct->l]);
  MPI qh = mpi_copy(hectx.qh[ct->l]);
  MPI Delta = mpi_new(0);
  s64_to_mpi(Delta, hectx.Delta);
  /* main */
  for (unsigned int i=0; i<n; i++) {
    mpi_rdiv(ct->c0.coeffs[i], ct->c0.coeffs[i], Delta);
    mpi_rdiv(ct->c1.coeffs[i], ct->c1.coeffs[i], Delta);
    mpi_smod(ct->c0.coeffs[i], q, qh);
    mpi_smod(ct->c1.coeffs[i], q, qh);
  }
  mpi_release(Delta);
  mpi_release(q);
  mpi_release(qh);
}

END_DECLS
