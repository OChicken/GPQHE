/*
 * Sample.
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
#include "poly.h"
#include <stddef.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>

BEGIN_DECLS

/* poly.h */
extern struct poly_ctx polyctx;

/* rng.c */
extern void randombytes_init(uint8_t *entropy_input, uint8_t *personalization_string);
extern void randombytes(uint8_t *x,size_t xlen);

void sample_z01vec(_Complex double vec[], const size_t m)
{
  uint8_t buf[m*2];
  randombytes(buf, 2*m);
  for (size_t i=0; i<m; i++)
    vec[i] = (double)buf[i]/256 + I*(double)buf[i+m]/256;
}

void sample_discrete_gaussian(int64_t *vec, const size_t m)
{
  uint8_t buf[m];
  randombytes(buf, m);
  for (size_t i=0; i<m; i+=2) {
    double r1 = (double)buf[i  ]/256;
    double r2 = (double)buf[i+1]/256;
    double theta = 2*GPQHE_PI*r1;
    double rr = sqrt(-2*log(r2)) * GPQHE_SIGMA;
    vec[i  ] = (int16_t)floor(rr*cos(theta)+0.5);
    vec[i+1] = (int16_t)floor(rr*sin(theta)+0.5);
  }
}

void sample_hwt(int64_t *vec, const size_t m)
{
  uint8_t tmp[GPQHE_BLKSIZ/8];
  randombytes(tmp, sizeof(tmp)); /* Generate noise in blocks of 64 coefficients */
  uint64_t num = loadu64_littleendian(tmp);
  unsigned int logm = 63-__builtin_clzll(m);
  size_t idx = 0; /* total nonzero entries count */
  while (idx<GPQHE_BLKSIZ) {
    uint8_t buf[8];
    randombytes(buf, 8);
    uint64_t i = loadnbits_littleendian(buf, logm);
    if (vec[i]==0) {
      vec[i] = (((num>>idx)&0x01)==0)? 1 : -1;
      idx++;
    }
  }
}

void sample_zero_center(int64_t *vec, const size_t m)
{
  uint8_t buf[2*m/8];
  randombytes(buf, sizeof(buf));
  MPI a = gcry_mpi_new(0);
  loadmpi_littleendian(a, buf, 2*m);
  for (size_t i=0; i<m; i++)
    vec[i] = (mpi_test_bit(a, 2*i)==0)? 0 : (mpi_test_bit(a, 2*i+1)==0)? 1: -1;
  mpi_release(a);
}

void sample_uniform(poly_mpi_t *r, const MPI q)
{
  unsigned int qbits = mpi_get_nbits(q);
  for (size_t i=0; i<polyctx.n; i++){
    uint8_t buf[qbits/8+1];
    randombytes(buf, sizeof(buf));
    loadmpi_littleendian(r->coeffs[i], buf, qbits);
  }
}

END_DECLS
