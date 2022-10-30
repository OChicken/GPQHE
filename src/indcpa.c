/*
 * indcpa - INDistinguishability under Chosen-Plaintext Attack.
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
#include "kem.h"

BEGIN_DECLS

/* poly.h */
extern struct poly_ctx polyctx;

/* kem.h */
extern struct kem_ctx kemctx;

/* rng.c */
extern void randombytes(uint8_t *x,size_t xlen);

#if 0
/**
 * indcpa_keypair - Generates public and private key for the CPA public-key
 * encryption scheme underlying the GPQHE.
 *
 * @pk[]: pointer to output public key
 *        (an already allocated array of CRYPTO_PUBLICKEYBYTES bytes)
 * @sk[]: pointer to output private key
 *        (an already allocated array of CRYPTO_SECRETKEYBYTES bytes)
 */
void indcpa_keypair(uint8_t pk[], uint8_t sk[])
{
  poly_mpi_t s, e;
  poly_rns_t shat, ehat;
  uint8_t z[2*GPQHE_SYMBYTES];
  uint8_t *publicseed = z;
  uint8_t nonce = 0;
  uint8_t *noiseseed = z+GPQHE_SYMBYTES;

  z[0] = 0x01;
  randombytes(z+1, GPQHE_SYMBYTES);
  shake256(z, 2*GPQHE_SYMBYTES, z, GPQHE_SYMBYTES + 1);

  poly_mpi_alloc(&s);
  poly_mpi_alloc(&e);

  poly_sample(&s, noiseseed, nonce++); /* nonce=0 */
  poly_sample(&e, noiseseed, nonce++); /* nonce=1 */

  poly_rns_alloc(&shat);
  poly_rns_alloc(&ehat);

  poly_rns_free(&shat);
  poly_rns_free(&ehat);

  poly_mpi_free(&s);
  poly_mpi_free(&e);
}


void indcpa_enc(uint8_t ct[],//GPQHE_INDCPA_BYTES
                const uint8_t m[],//GPQHE_INDCPA_MSGBYTES
                const uint8_t pk[],//GPQHE_INDCPA_PUBLICKEYBYTES
                const uint8_t coins[GPQHE_SYMBYTES])
{
  //poly_rns_t sp;
  uint8_t nonce = 0;
  poly_mpi_t sp, ep, epp;

  poly_mpi_alloc(&sp);
  poly_mpi_alloc(&ep);
  poly_mpi_alloc(&epp);

  poly_sample(&sp,  coins, nonce++); /* nonce=0 */
  poly_sample(&ep,  coins, nonce++); /* nonce=1 */
  poly_sample(&epp, coins, nonce++); /* nonce=2 */

  poly_mpi_free(&sp);
  poly_mpi_free(&ep);
  poly_mpi_free(&epp);
}
#endif

END_DECLS
