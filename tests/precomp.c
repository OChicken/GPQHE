/*
 * Test file: precompute constants and check with HEAAN.
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

static void test_precomp()
{
TEST_BEGIN();

  unsigned int logq = 800;
  MPI q = mpi_new(0);
  mpi_set_ui(q, 1);
  mpi_lshift(q, q, logq);
TEST_DO("init");
  polyctx_init(16,q);
TEST_DONE();

  uint64_t ptr;

TEST_DO("check ntt parameters");
  int fd_pVec, fd_pInvVec, fd_prVec;
  int fd_ninv, fd_zetas, fd_zetas_inv;
  fd_pVec      = open("HEAAN_KAT/pVec",      O_RDONLY, 0000);
  fd_pInvVec   = open("HEAAN_KAT/pInvVec",   O_RDONLY, 0000);
  fd_prVec     = open("HEAAN_KAT/prVec",     O_RDONLY, 0000);
  fd_ninv      = open("HEAAN_KAT/ninv",      O_RDONLY, 0000);
  fd_zetas     = open("HEAAN_KAT/zetas",     O_RDONLY, 0000);
  fd_zetas_inv = open("HEAAN_KAT/zetas_inv", O_RDONLY, 0000);
  for (struct rns_ctx *rns = polyctx.rns; rns; rns=rns->next) {
    read(fd_pVec, &ptr, sizeof(uint64_t));
    assert(rns->p==ptr);
    read(fd_pInvVec, &ptr, sizeof(uint64_t));
    assert(rns->pinv_mont==ptr);
    read(fd_prVec, &ptr, sizeof(uint64_t));
    assert(rns->pinv_barr==ptr);
    read(fd_ninv, &ptr, sizeof(uint64_t));
    assert(rns->ninv==ptr);
    for (unsigned int i=0; i<polyctx.n; i++) {
      read(fd_zetas, &ptr, sizeof(uint64_t));
      assert(rns->zetas[i]==ptr);
      read(fd_zetas_inv, &ptr, sizeof(uint64_t));
      assert(rns->zetas_inv[i]==ptr);
    }
  }
  close(fd_pVec);
  close(fd_pInvVec);
  close(fd_prVec);
  close(fd_ninv);
  close(fd_zetas);
  close(fd_zetas_inv);
TEST_DONE();

TEST_DO("check crt parameters");
  int fd_phat_invmp;
  fd_phat_invmp = open("HEAAN_KAT/phat_invmp", O_RDONLY, 0000);
  
  for (struct rns_ctx *rns=polyctx.rns; rns; rns=rns->next) {
    for (unsigned int d=0; d<rns->dim; d++) {
      read(fd_phat_invmp, &ptr, sizeof(uint64_t));
      assert(rns->phat_invmp[d]==ptr);
    }
  }
  close(fd_phat_invmp);
TEST_DONE();

TEST_DO("exit");
  polyctx_exit();
TEST_DONE();

  gcry_mpi_release(q);
TEST_END();
}

int main()
{
  test_precomp();
  return 0;
}
