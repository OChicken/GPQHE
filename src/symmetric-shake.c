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

#include <stddef.h>
#include <stdint.h>
#include "params.h"
#include "fips202.h"
#include "symmetric.h"

BEGIN_DECLS

void hash_h(uint8_t h[32], const uint8_t *in, size_t inlen)
{
  sha3_256(h, in, inlen);
}

void hash_g(uint8_t h[64], const uint8_t *in, size_t inlen)
{
  sha3_512(h, in, inlen);
}

/**
 * xof_absorb - Absorb step of the SHAKE128 specialized for the GPQHE context.
 * 
 * @state: pointer to (uninitialized) output Keccak state
 * @seed: pointer to GPQHE_SYMBYTES input to be absorbed into state
 * @x: additional byte of input
 * @y: additional byte of input
 */
void xof_absorb(keccak_state *state,
                const uint8_t seed[GPQHE_SYMBYTES],
                uint8_t x, uint8_t y)
{
  uint8_t extseed[GPQHE_SYMBYTES+2];

  memcpy(extseed, seed, GPQHE_SYMBYTES);
  extseed[GPQHE_SYMBYTES+0] = x;
  extseed[GPQHE_SYMBYTES+1] = y;

  shake128_absorb_once(state, extseed, sizeof(extseed));
}

void xof_squeezeblocks(uint8_t *out, size_t nblocks, keccak_state *state)
{
  shake128_squeezeblocks(out, nblocks, state);
}

/**
 * gpqhe_shake256_prf - Usage of SHAKE256 as a PRF, concatenates secret and
 * public input and then generates outlen bytes of SHAKE256 output
 *
 * @out: pointer to output
 * @outlen: number of requested output bytes
 * @key: pointer to the key (of length GPQHE_SYMBYTES)
 * @nonce: single-byte nonce (public PRF input)
 */
void prf(uint8_t *out, size_t outlen, const uint8_t key[GPQHE_SYMBYTES], uint8_t nonce)
{
  uint8_t extkey[GPQHE_SYMBYTES+1];

  memcpy(extkey, key, GPQHE_SYMBYTES);
  extkey[GPQHE_SYMBYTES] = nonce;

  shake256(out, outlen, extkey, sizeof(extkey));
}

void kdf(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen)
{
  shake256(out, outlen, in, inlen);
}

END_DECLS
