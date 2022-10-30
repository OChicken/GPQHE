/*
 * FHE operation collections of Galois automorphism: conjugate and rotate.
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
extern void mpi_smod(MPI r, const MPI q, const MPI qh);

/** Key-switching procedure */
static void he_swk(he_ct_t *ct,
  const poly_mpi_t *d0, const poly_mpi_t *d1,
  const he_pk_t *swk, const MPI ql)
{
  /* local variables */
  MPI qlh = mpi_copy(hectx.qh[ct->l]);
  MPI P   = mpi_copy(hectx.P);
  MPI Pql = mpi_new(0);
  mpi_mul(Pql, P, ql);
  unsigned int n = polyctx.n;
  unsigned int dim = (mpi_get_nbits(ql)+mpi_get_nbits(Pql))/GPQHE_LOGP+1;
  /* main */
  poly_mul(&ct->c0, d1, &swk->p0, dim, Pql); /* d0 in Rql, swk.p0 in R_Pql */
  poly_mul(&ct->c1, d1, &swk->p1, dim, Pql); /* d1 in Rql, swk.p1 in R_Pql */
  for (unsigned int i=0; i<n; i++) {
    mpi_rdiv(ct->c0.coeffs[i], ct->c0.coeffs[i], P); /* (d1*swk.p0)/P, now in Rql */
    mpi_rdiv(ct->c1.coeffs[i], ct->c1.coeffs[i], P); /* (d1*swk.p1)/P, now in Rql */
    mpi_addm(ct->c0.coeffs[i], ct->c0.coeffs[i], d0->coeffs[i], ql);
    mpi_smod(ct->c0.coeffs[i], ql, qlh);
    mpi_smod(ct->c1.coeffs[i], ql, qlh);
  }
  /* release */
  mpi_release(qlh);
  mpi_release(P);
  mpi_release(Pql);
}

void he_conj(he_ct_t *ct, const he_pk_t *ck)
{
  /* local varables */
  MPI q = mpi_copy(hectx.q[ct->l]);
  poly_mpi_t d0, d1;
  poly_mpi_alloc(&d0);
  poly_mpi_alloc(&d1);
  /* main */
  poly_conj(&d0, &ct->c0);
  poly_conj(&d1, &ct->c1);
  he_swk(ct, &d0, &d1, ck, q);
  /* release */
  mpi_release(q);
  poly_mpi_free(&d0);
  poly_mpi_free(&d1);
}

void he_rot(he_ct_t *ct, const int rot, const he_pk_t *rk)
{
  /* local variables */
  MPI q = mpi_copy(hectx.q[ct->l]);
  poly_mpi_t d0, d1;
  poly_mpi_alloc(&d0);
  poly_mpi_alloc(&d1);
  /* main */
  poly_rot(&d0, &ct->c0, rot);
  poly_rot(&d1, &ct->c1, rot);
  he_swk(ct, &d0, &d1, &rk[rot], q);
  /* release */
  mpi_release(q);
  poly_mpi_free(&d0);
  poly_mpi_free(&d1);
}

END_DECLS
