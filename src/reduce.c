/*
 * Montgomery reduction and Barrett reduction.
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
#include <errno.h>

BEGIN_DECLS

/* poly.h */
extern struct poly_ctx polyctx;

/**
 * montgomery_inv - get the inverse of input q in the Montgomery domain
 * 
 * @q: the 64-bit modulue
 * 
 * return qinv, the inverse of q in the Montgomery domain 2^64
 */
uint64_t montgomery_inv(uint64_t q)
{
  int64_t qinv = 1;
  for (unsigned int k=64; k>0; k--) {
    qinv *= q;
    qinv &= polyctx.Rsub1;
    q    *= q;
    q    &= polyctx.Rsub1;
  }
  //if (qinv >= (1ULL<<32)) /* qinv >= R/2 */
  //  qinv = (int64_t)qinv; /* qinv -= R */
  return qinv;
}

/**
 * montgomery_reduce - Multiplication followed by Montgomery reduction
 * 
 * @a: first factor
 * @q: the modulue
 * @qinv: the inverse of q in the Montgomery domain
 * 
 * return a 64-bit integer congruent to a*R^{-1} mod q
 */
uint64_t montgomery_reduce(u128 a, uint64_t q, int64_t qinv)
{
  uint64_t lo =  a               & polyctx.Rsub1; /* a mod R, >0 */
  uint64_t hi = (a>>polyctx.logR)& polyctx.Rsub1;
  uint64_t u  = ((u128)lo * (u128)qinv) & polyctx.Rsub1;
  uint64_t t  = ((u128)u  * (u128)q)   >> polyctx.logR;
  return (hi<t) ? hi-t+q : hi-t;
}

/**
 * barrett_inv - get the inverse of input q in the Barrett domain
 * 
 * @q: the 64-bit modulue
 * 
 * return qinv = 2^(2*nbits(q)) / q
 */
uint64_t barrett_inv(uint64_t q)
{
  return ((u128)1 << (2*(64-__builtin_clzll(q)))) / (u128)q;
}

/**
 * barrett_reduce - get the inverse of input q in the Barrett domain
 * 
 * @q: the 64-bit modulue
 * @qinv: 2^(2*nbits(q)) / q
 * 
 * return a mod q
 */
uint64_t barrett_reduce(u128 a, uint64_t q, uint64_t qinv)
{
  uint64_t lo =  a                & polyctx.Rsub1; /* a mod R, >0 */
  uint64_t hi = (a>>polyctx.logR) & polyctx.Rsub1;
  u128 t = (((u128)lo * (u128)qinv)>>polyctx.logR)
           +((u128)hi * (u128)qinv);
  int shift = 2*(64-__builtin_clzll(q)) - 64;
  if (shift<0)  {
    errno = EINVAL;
    fprintf(stderr, "\033[1m\033[31merror:\033[0m \033[1m%s\033[0m. "
      "The number of bits of the modulus is too small.\n", strerror(errno));
    abort();
  }
  assert(shift>=0);
  t >>= shift;
  a -= t*q;
  uint64_t r = a & polyctx.Rsub1;
  return (r<q)? r : r-q;
}

END_DECLS
