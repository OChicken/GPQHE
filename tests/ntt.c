/*
 * Test file: ntt check with HEAAN.
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
#include <assert.h>
#include <fcntl.h>   /* open & creat, return int */
#include <unistd.h>  /* read & write, return ssize_t */
#include <pmu.h>

/* poly.h */
extern struct poly_ctx polyctx;

/* ntt.c */
extern void    ntt(uint64_t a[], const struct rns_ctx *rns);
extern void invntt(uint64_t a[], const struct rns_ctx *rns);

static void test_ntt()
{
TEST_BEGIN();

  unsigned int logq = 800;
  MPI q = mpi_new(0);
  mpi_set_ui(q, 1);
  mpi_lshift(q, q, logq);
TEST_DO("init");
  polyctx_init(16,q);
TEST_DONE();

  uint64_t ptr[polyctx.n];
  int fd_ntt = open("HEAAN_KAT/ntt", O_RDONLY, 0000);

  uint64_t rai0[polyctx.n];
  uint64_t rai[polyctx.n];
  for (unsigned int i=0; i<polyctx.n; i++)
    rai0[i] = rai[i] = i;

TEST_DO("ntt");
  
  for (struct rns_ctx *rns = polyctx.rns; rns; rns=rns->next) {
    read(fd_ntt, ptr, polyctx.n*sizeof(uint64_t));
    ntt(rai, rns);
    for (unsigned int i=0; i<polyctx.n; i++)
      assert(rai[i]==ptr[i]);
    invntt(rai, rns);
    for (unsigned int i=0; i<polyctx.n; i++)
      assert(rai[i]==rai0[i]);
  }
TEST_DONE();
  close(fd_ntt);

TEST_DO("exit");
  polyctx_exit();
TEST_DONE();

  gcry_mpi_release(q);
TEST_END();
}

int main()
{
  test_ntt();
  return 0;
}
