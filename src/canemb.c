/*
 * Canonical embedding.
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

/* polyctx.h */
extern struct poly_ctx polyctx;

static void bitrev_vec(_Complex double a[], const unsigned int n)
{
  for (unsigned int i=1,j=0; i<n; i++) {
    unsigned int bit = n>>1;
    for (; j>=bit; bit>>=1)
      j -= bit;
    j += bit;
    if (i<j) {
      _Complex double tmp = a[i];
      a[i]=a[j];
      a[j]=tmp;
    }
  }
}

void canemb(_Complex double a[], const unsigned int slots)
{
  bitrev_vec(a, slots);
  for (unsigned int len=2; len<=slots; len<<=1) {
    unsigned int idx_mod = len<<2;
    unsigned int gap = polyctx.m/idx_mod;
    unsigned int mid = len>>1;
    for (unsigned int i=0; i<slots; i+=len) {
      for (unsigned int j=0; j<mid; j++) {
        unsigned int k = (polyctx.ring.cyc_group[j] % idx_mod) * gap;
        _Complex double u = a[i+j];
        _Complex double v = a[i+j+mid]*polyctx.ring.zetas[k];
        a[i+j]     = u+v;
        a[i+j+mid] = u-v;
      }
    }
  }
}

void invcanemb(_Complex double a[], const unsigned int slots)
{
  for (unsigned int len=slots; len>=1; len>>=1) {
    unsigned int idx_mod = len << 2;
    unsigned int gap = polyctx.m/idx_mod;
    unsigned int mid = len >> 1;
    for (unsigned int i=0; i<slots; i+=len) {
      for (unsigned int j=0; j<mid; j++) {
        unsigned int k = (idx_mod - (polyctx.ring.cyc_group[j] % idx_mod)) * gap;
        _Complex double u =  a[i+j] + a[i+j+mid];
        _Complex double v = (a[i+j] - a[i+j+mid])*polyctx.ring.zetas[k];
        a[i+j]     = u;
        a[i+j+mid] = v;
      }
    }
  }
  bitrev_vec(a, slots);
  for (unsigned int i=0; i<slots; i++)
    a[i] /= slots;
}

END_DECLS
