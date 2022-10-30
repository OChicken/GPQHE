/*
 * Polynomial ring definitions and parameters.
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

#ifndef POLY_H
#define POLY_H

#include "params.h"

BEGIN_DECLS

struct rns_ctx {
  unsigned int dim; /* current number of limbs; current dimension of the RNS */
  uint64_t p;    /* prime base; prime limb; prime modulus */
  uint64_t pinv_mont;
  uint64_t pinv_barr;
  uint64_t ninv;
  uint64_t *zetas;
  uint64_t *zetas_inv;
  MPI P;      /* base product P */
  MPI P_2;    /* P/2 */
  MPI *phat;  /* phat_i = (P/pi) mod pi */
  uint64_t *phat_invmp; /* invm(phat%p, p) = phat^{-1} mod p */
  struct rns_ctx *next;
};

struct ring_ctx {
  unsigned int *cyc_group;
  _Complex double *zetas;
  double cM; /* ring constant */
};

struct poly_ctx {
  /* polynomial degree */
  unsigned int logn;
  unsigned int n;
  unsigned int m;
  /* modulus q */
  unsigned int logqub; /* upper bound of logq, depends on logn by table look-up */
  unsigned int logq;   /* log(q) */
  MPI q; /* modulus must guarantee logq <= logqub ; current modulus */
  /* Montgomery R */
  unsigned int logR;
  u128 R;
  u128 Rsub1;
  unsigned int dimub; /* maximum number of limbs; maximum dimension of the RNS */
  struct rns_ctx *rns;   /* RNS parameters */
  struct ring_ctx ring; /* ring parameters */
};

/* Elements of R_q = Z_q[X]/(X^n + 1). Represents polynomial
 * coeffs[0] + X*coeffs[1] + X^2*xoeffs[2] + ... + X^{n-1}*coeffs[n-1] */
struct poly_mpi {
  MPI *coeffs;
};
typedef struct poly_mpi poly_mpi_t;

struct poly_rns {
  uint64_t *coeffs;
};
typedef struct poly_rns poly_rns_t;

/* poly.c */
void poly_mpi_alloc(poly_mpi_t *a);
void poly_mpi_free(poly_mpi_t *a);
void poly_rns_alloc(poly_rns_t *a, const unsigned int dim);
void poly_rns_free(poly_rns_t *a);
void poly_rns_add(uint64_t rhat[], const uint64_t ahat[], const uint64_t bhat[], const struct rns_ctx *rns);
void poly_rns_mul(uint64_t rhat[], const uint64_t ahat[], const uint64_t bhat[], const struct rns_ctx *rns);
void poly_mul(poly_mpi_t *r, const poly_mpi_t *a, const poly_mpi_t *b,
              const unsigned int dim, const MPI q);
void poly_rns2mpi(poly_mpi_t *r, const poly_rns_t *rhat, const struct rns_ctx *rns, const MPI q);
void poly_rot(poly_mpi_t *r, const poly_mpi_t *a, size_t rot);
void poly_conj(poly_mpi_t *r, const poly_mpi_t *a);
void poly_sample(poly_mpi_t *r, const unsigned char *seed, unsigned char nonce);

/* precomp.c */
void polyctx_init(unsigned int logn, MPI q);
void polyctx_exit();

END_DECLS

#endif /* POLY_H */
