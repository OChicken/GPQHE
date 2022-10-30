/*
 * NTT number theoretic transform.
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

BEGIN_DECLS

/* poly.h */
extern struct poly_ctx polyctx;

/* reduce.c */
extern uint64_t montgomery_reduce(u128 a, uint64_t q, int64_t qinv);
extern uint64_t barrett_reduce(u128 a, uint64_t q, uint64_t qinv);

static uint64_t fqmul(uint64_t a, uint64_t b, uint64_t q, uint64_t qinv)
{
  return montgomery_reduce((u128)a*(u128)b, q, qinv);
}

void ntt(uint64_t a[], const struct rns_ctx *rns)
{
  const uint64_t q    = rns->p;
  const uint64_t qinv = rns->pinv_mont;
  unsigned int len, start, j=0, k=1;
  for (len=polyctx.n/2; len>=1; len>>=1) { /* len=N/2,...,1 */
    for (start=0; start<polyctx.n; start=j+len) {
      uint64_t zeta = rns->zetas[k++];
      for (j=start; j<start+len; j++) {
        uint64_t t = fqmul(a[j+len], zeta, q, qinv);
        a[j+len] = (a[j]>=t  )? a[j]-t : a[j]-t+q;
        a[j]     = (a[j]<=q-t)? a[j]+t : a[j]+t-q;
      }
    }
  }
}

void invntt(uint64_t a[], const struct rns_ctx *rns)
{
  const uint64_t q    = rns->p;
  const uint64_t qinv = rns->pinv_mont;
  const uint64_t ninv = rns->ninv;
  unsigned int len, start, j=0, k=0;
  for (len=1, k=polyctx.n/(2*len); len<=polyctx.n/2; len<<=1, k=polyctx.n/(2*len)) {
    for (start=0; start<polyctx.n; start=j+len) {
      uint64_t zeta_inv = rns->zetas_inv[k++];
      for (j=start; j<start+len; j++) {
        uint64_t t = a[j];
        a[j]     = (a[j+len]<=q-t)? t+a[j+len] : t+a[j+len]-q;
        a[j+len] = (a[j+len]<=t  )? t-a[j+len] : t-a[j+len]+q;
        a[j+len] = fqmul(a[j+len], zeta_inv, q, qinv);
      }
    }
  }
  for (long i = 0; i < polyctx.n; i++)
    a[i] = fqmul(a[i], ninv, q, qinv);
}

END_DECLS
