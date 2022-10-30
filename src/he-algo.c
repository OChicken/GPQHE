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

static void zrotdiag(_Complex double rotdiag[], const _Complex double A[],
                          const unsigned int idx, const int rot)
{
  _Complex double diag[hectx.slots];
  const unsigned int m = hectx.slots;
  for (unsigned int i=0; i<m; i++)
    diag[i] = A[(i%m)*m + (idx+i)%m];
  for (unsigned int i=0; i<m; i++) {
    int idx = ((int)i+rot)%m;
    if (idx<0)
      idx += m;
    assert(idx>=0);
    rotdiag[i] = diag[idx];
  }
}

/** gemv: GEneral Matrix Vector multiplication
 * he_ckks_gemv(&ct_Av, A, &ct_v, evk->rot_key, ctx); */
void he_gemv(he_ct_t *ct_dest, _Complex double *A, const he_ct_t *ct, const he_pk_t *rk)
{
  /* plaintext A is of dim (d/2)*(d/2).  */
  unsigned int n1 = sqrt(hectx.slots);
  if (hectx.slots!=n1*n1)
    n1 = sqrt(2*hectx.slots);
  unsigned int n2 = hectx.slots/n1;

  /* alloc */
  he_pt_t pt;
  he_alloc_pt(&pt);
  he_ct_t inner_sum, ct_rot;
  he_alloc_ct(&inner_sum);
  he_alloc_ct(&ct_rot);
  /* Compute sum */
  for (unsigned int i=0; i<n2; i++) {
    int shift = i*n1;
    for (unsigned int j=0; j<n1; j++) {
      /* rotate and encode ct */
      he_copy_ct(&ct_rot, ct);
      he_rot(&ct_rot, j, rk);
      /* rotate and encode A */
      _Complex double rotdiag[hectx.slots];
      zrotdiag(rotdiag, A, shift+j, -shift);
      he_ecd(&pt, rotdiag);
      /* multiply and accumulate */
      he_mulpt(&ct_rot, &ct_rot, &pt);
      if (!j)
        he_copy_ct(&inner_sum, &ct_rot);
      else
        he_add(&inner_sum, &inner_sum, &ct_rot);
    }
    he_rot(&inner_sum, shift, rk);
    if (!i)
      he_copy_ct(ct_dest, &inner_sum);
    else
      he_add(ct_dest, ct_dest, &inner_sum);
  }
  he_rs(ct_dest);
  /* release */
  he_free_pt(&pt);
  he_free_ct(&inner_sum);
  he_free_ct(&ct_rot);
}


END_DECLS
