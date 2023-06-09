/*
 * Symmetric encrytion.
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

#include "params.h"
#include "fips202.h"

#ifndef SYMMETRIC_H
#define SYMMETRIC_H

BEGIN_DECLS

typedef keccak_state xof_state;

#define XOF_BLOCKBYTES SHAKE128_RATE

void hash_h(uint8_t h[32], const uint8_t *in, size_t inlen);
void hash_g(uint8_t h[64], const uint8_t *in, size_t inlen);
#define xof_absorb gpqhe_shake128_xof_absorb
void xof_absorb(keccak_state *state, const uint8_t seed[GPQHE_SYMBYTES], uint8_t x, uint8_t y);
void xof_squeezeblocks(uint8_t *out, size_t nblocks, keccak_state *state);
#define prf gpqhe_shake256_prf
void prf(uint8_t *out, size_t outlen, const uint8_t key[GPQHE_SYMBYTES], uint8_t nonce);
#define kdf gpqhe_shake256_kdf
void kdf(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen);

END_DECLS

#endif /* SYMMETRIC_H */
