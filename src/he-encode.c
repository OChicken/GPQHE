/*
 * Encoder.
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

#include "poly.h"
#include "fhe.h"
#include <complex.h> /* creal, cimag, I */
#include <math.h>    /* round */
#include <assert.h>

BEGIN_DECLS

/* polyctx.h */
extern struct poly_ctx polyctx;

/* fhe.h */
extern struct he_ctx hectx;

/* canemb.c */
void canemb(_Complex double a[], const unsigned int slots);
void invcanemb(_Complex double a[], const unsigned int slots);

/** Encodes complex array into pt, an integral polynomial. */
void encode(poly_mpi_t *r, const _Complex double m[], const double Delta)
{
  _Complex double u[hectx.slots];
  memcpy(u, m, sizeof(u));
  unsigned int nh = polyctx.n/2; /* n half */
  unsigned int gap = nh/hectx.slots;
  invcanemb(u, hectx.slots);
  for (unsigned int i=0, j=0; i<hectx.slots; i++, j+=gap) {
    double_to_mpi(&r->coeffs[j   ], round((long double)creal(u[i])*Delta));
    double_to_mpi(&r->coeffs[j+nh], round((long double)cimag(u[i])*Delta));
  }
}

/** Decodes pt, an integral polynomial, into complex array. */
void decode(_Complex double m[], poly_mpi_t *a, const double Delta, const MPI q)
{
  MPI real = mpi_new(0);
  MPI imag = mpi_new(0);
  unsigned int nh = polyctx.n/2;
  unsigned int gap = nh/hectx.slots;
  unsigned int logq = mpi_get_nbits(q);
  for (unsigned int i=0, j=0; i<hectx.slots; i++, j+=gap) {
    mpi_set(real, a->coeffs[j   ]);
    mpi_set(imag, a->coeffs[j+nh]);
#if 0
    mpi_mod(real, a->coeffs[j   ], q);
    mpi_mod(imag, a->coeffs[j+nh], q);
    if (mpi_get_nbits(real)==logq)
      mpi_sub(real, real, q);
    if (mpi_get_nbits(imag)==logq)
      mpi_sub(imag, imag, q);
#endif
    m[i] = (double)mpi_to_s128(real)/Delta
        +I*(double)mpi_to_s128(imag)/Delta;
  }
  canemb(m, hectx.slots);
  mpi_release(real);
  mpi_release(imag);
}

/** canemb_norm - Canonical embedding norm */
double canemb_norm(const _Complex double m[], const double Delta)
{
  unsigned int len = hectx.slots*2;
  double u[len];
  for (unsigned int i=0; i<hectx.slots; i++) {
    u[i            ] = round((long double)creal(m[i])*Delta);
    u[i+hectx.slots] = round((long double)cimag(m[i])*Delta);
  }
  double norm = 0;
  for (unsigned int i=0; i<len; i++) {
    double val = fabs(u[i]);
    if (norm<val)
      norm = val;
  }
  return norm;
}

/** Encodes complex array into pt, an integral polynomial. */
void he_ecd(struct he_pt *pt, const _Complex double *m)
{
  pt->l = hectx.L;
  double norm = canemb_norm(m, hectx.Delta);
  pt->nu = (norm>=hectx.Delta)? norm : hectx.Delta;
  encode(&pt->m, m, pt->nu);
}

/** Decodes pt, an integral polynomial, into complex array. */
void he_dcd(_Complex double *m, struct he_pt *pt)
{
  double Delta = (pt->nu>=hectx.Delta)? hectx.Delta : pt->nu;
  decode(m, &pt->m, Delta, hectx.q[pt->l]);
}

END_DECLS
