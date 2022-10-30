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

#ifndef KEM_H
#define KEM_H

#include "poly.h"

BEGIN_DECLS

struct kem_ctx {
  unsigned int polybytes; /* bytes required to represent a polynomial */
  unsigned int pkbytes;   /* public key bytes */
  unsigned int skbytes;   /* secret key bytes */
  unsigned int ctbytes;   /* ciphertext bytes*/
  unsigned int ssbytes;   /* size in bytes of shared key */
};

/* precomp.c */
void kemctx_init(unsigned int ssbytes);
void kemctx_exit();

END_DECLS

#endif /* KEM_H */
