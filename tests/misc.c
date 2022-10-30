/*
 * Test file: misc.
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
#include "fhe.h"
#include <assert.h>
#include <fcntl.h>   /* open & creat, return int */
#include <unistd.h>  /* read & write, return ssize_t */
#include <pmu.h>

#include <complex.h>
#include <math.h>

/* poly.h */
extern struct poly_ctx polyctx;

/* fhe.h */
extern struct he_ctx hectx;

/* sample.c */
extern void sample_z01vec(_Complex double vec[], const unsigned int m);
extern void sample_hwt(int64_t *vec, const size_t m);
extern void sample_zero_center(int64_t *vec, const size_t m);

static void test_misc()
{
TEST_BEGIN();
  MPI q=mpi_new(0);
  mpi_set_ui(q, 1);
  mpi_lshift(q, q, 220);
  unsigned int logn = 14;

TEST_DO("init");
  polyctx_init(logn, q);
TEST_DONE();

  uint8_t *buf;
  buf=malloc(2);
  buf[0]=0x32;
  buf[1]=0x54;
  printf("%lu\n", loadnbits_littleendian(buf,11));
  loadmpi_littleendian(q, buf, 11);
  show_mpi(q);
  printf("%lu\n", mpi_to_u64(q));
  free(buf);

  buf=malloc(9);
  for (unsigned int i=0; i<8; i++)
    buf[i]=0xff;
  buf[8]=0x55;
  loadmpi_littleendian(q, buf, 66);
  show_mpi(q);
  free(buf);

TEST_DO("loadu64/uint_littleendian");{
  uint8_t *buf;
  buf=malloc(8);
  buf[0]=0xc6;
  buf[1]=0x23;
  buf[2]=0x7b;
  buf[3]=0x32;
  buf[4]=0x67;
  buf[5]=0x45;
  buf[6]=0x8b;
  buf[7]=0x6b;
  uint32_t hi = 1804289383;//rand() first call
  uint32_t lo = 846930886; //rand() second call
  uint64_t comb = ((uint64_t)hi<<32)+lo;
  TEST_PLACEHOLDER();
  printf("hi=%u lo=%u\n", hi, lo);
  TEST_PLACEHOLDER();
  printf("hi<<32+lo = %lu = %#02lx\n", comb, comb);
  assert(loadu64_littleendian(buf)==comb);
  assert(loadnbits_littleendian(buf,64)==comb);
  free(buf);
}TEST_DONE();

  mpi_release(q);
  polyctx_exit();
TEST_END();
}

__attribute__((unused))
static void test_sample()
{
  unsigned int len=16;
  _Complex double zvec[len];
  sample_z01vec(zvec, len);
  for (unsigned int i=0; i<len; i++)
    printf("(%f,%f)\n", creal(zvec[i]), cimag(zvec[i]));
  unsigned int zeros, pones, nones;
  len = 65536;
  int64_t vec[len];
  memset(vec, 0, sizeof(vec));
  sample_hwt(vec, len);
  zeros = pones = nones = 0;
  for (unsigned int i=0; i<len; i++) {
    if (vec[i]== 0) zeros++;
    if (vec[i]==+1) pones++;
    if (vec[i]==-1) nones++;
  }
  printf("0's=%u, +1's=%u, -1's=%u\n", zeros, pones, nones);
  memset(vec, 0, sizeof(vec));
  sample_zero_center(vec, len);
  zeros = pones = nones = 0;
  for (unsigned int i=0; i<len; i++) {
    if (vec[i]== 0) zeros++;
    if (vec[i]==+1) pones++;
    if (vec[i]==-1) nones++;
  }
  printf("0's=%u, +1's=%u, -1's=%u\n", zeros, pones, nones);
}

int main()
{
  test_misc();
  //test_sample();
  return 0;
}
