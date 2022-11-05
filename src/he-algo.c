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

/* fhe.h */
extern struct he_ctx hectx;

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
void he_gemv(he_ct_t *ct_dest, _Complex double *A, const he_ct_t *ct, const he_evk_t *rk)
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

void he_inv(he_ct_t *ct_inv, const he_ct_t *ct, const he_evk_t *rlk, const unsigned int iter)
{
  /* alloc */
  he_pt_t one, two;
  he_alloc_pt(&one);
  he_alloc_pt(&two);
  he_ct_t tmp, an, bn;
  he_alloc_ct(&tmp);
  he_alloc_ct(&an);
  he_alloc_ct(&bn);
  /* assign */
  he_const_pt(&one, 1);
  he_const_pt(&two, 2);
  /* main */
  he_copy_ct(&tmp, ct);
  he_neg(&tmp);              /* tmp = -ct */
  he_addpt(&an, &tmp, &two); /* an = 2-ct */
  he_moddown(&an);
  he_addpt(&bn, &tmp, &one); /* bn = 1-ct */
  for (unsigned int n=0; n<iter; n++) {
    he_mul(&bn, &bn, &bn, rlk);  /* bn = bn*bn */
    he_rs(&bn);
    he_addpt(&tmp, &bn, &one);         /* tmp = 1+bn */
    he_mul(&an, &an, &tmp, rlk);
    he_rs(&an);
  }
  he_copy_ct(ct_inv, &an);
  he_free_pt(&one);
  he_free_pt(&two);
  he_free_ct(&an);
  he_free_ct(&bn);
  he_free_ct(&tmp);
}

void he_sqrt(he_ct_t *ct_sqrt, const he_ct_t *ct, const he_evk_t *rlk, const unsigned int iter)
{
  /* alloc */
  he_pt_t one, three, half, quarter;
  he_alloc_pt(&one);
  he_alloc_pt(&three);
  he_alloc_pt(&half);
  he_alloc_pt(&quarter);
  he_ct_t tmp, an, bn;
  he_alloc_ct(&tmp);
  he_alloc_ct(&an);
  he_alloc_ct(&bn);
  /* assign */
  he_const_pt(&one, 1);
  he_const_pt(&three, 3);
  he_const_pt(&half, 0.5);
  he_const_pt(&quarter, 0.25);
  /* main */
  he_copy_ct(&an, ct);    /* an = ct */
  he_subpt(&bn, ct, &one); /* bn = ct-1 */
  for (size_t n=0; n<iter; n++) {
    /* an */
    he_mulpt(&tmp, &bn, &half);  /* bn/2 */
    he_rs(&tmp);
    he_subpt(&tmp, &tmp, &one);  /* bn/2-1 */
    he_neg(&tmp);                /* tmp = 1-bn/2 */
    he_moddown(&an);
    he_mul(&an, &an, &tmp, rlk); /* an = an*(1-bn/2) */
    he_rs(&an);
    /* bn */
    he_subpt(&tmp, &bn, &three);
    he_mulpt(&tmp, &tmp, &quarter); /* tmp = (bn-3)/4 */
    he_rs(&tmp);
    he_mul(&bn, &bn, &bn, rlk);     /* bn = bn*bn */
    he_rs(&bn);
    he_mul(&bn, &bn, &tmp, rlk);    /* bn = (bn*bn)*((bn-3)/4) */
    he_rs(&bn);
    he_show_ct_params(&an, "an[%u]", n);
  }
  he_copy_ct(ct_sqrt, &an);
}

END_DECLS
