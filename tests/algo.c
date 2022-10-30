/*
 * Test file: FHE algorithms collections.
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
#include <pmu.h>
#include <complex.h>
#include <math.h>

/* poly.h */
extern struct poly_ctx polyctx;

/* fhe.h */
extern struct he_ctx hectx;

/* sample.c */
extern void sample_z01vec(_Complex double vec[], const unsigned int m);

static double dznrmmax_dist(_Complex double y[], _Complex double x[],
                            const unsigned int len)
{
  double tmp = 0;
  double dist = 0;
  for (unsigned int i=0; i<len; i++) {
    tmp = cabs(y[i] - x[i]);
    if (dist<tmp)
      dist = tmp;
  }
  return dist;
}

static void test_gemv(const unsigned int logn,
                      const MPI q,
                      const unsigned int slots,
                      const uint64_t Delta,
                      const char *key)
{
TEST_BEGIN();

TEST_DO("init");
  hectx_init(logn, q, slots, Delta);
TEST_DONE();

  he_pt_t pt;
  he_alloc_pt(&pt);
  he_ct_t ct_Av, ct_v;
  he_alloc_ct(&ct_Av);
  he_alloc_ct(&ct_v);
  he_pk_t pk, rk[hectx.slots];
  he_alloc_pk(&pk);
  for (unsigned int rot=0; rot<hectx.slots; rot++)
    he_alloc_pk(&rk[rot]);
  poly_mpi_t sk;
  he_alloc_sk(&sk);
  _Complex double m[slots], Av[slots], A[slots*slots], v[slots];
  sample_z01vec(v, slots);
  sample_z01vec(A, slots*slots);
  memset(Av, 0, sizeof(Av));
  for (unsigned int i=0; i<slots; i++)
    for (unsigned int j=0; j<slots; j++)
      Av[i] += A[i*slots+j]*v[j];
  he_keypair(&pk, &sk);
TEST_DO("gen rk");
  he_genrk(rk, &sk);
TEST_DONE();

  he_ecd(&pt, v);
  if (!strcmp(key, "sk"))
    he_enc_sk(&ct_v, &pt, &sk);
  if (!strcmp(key, "pk"))
    he_enc_pk(&ct_v, &pt, &pk);
TEST_DO("gemv");
  he_gemv(&ct_Av, A, &ct_v, rk);
TEST_DONE();
  he_dec(&pt, &ct_Av, &sk);
  he_dcd(m, &pt);
  TEST_PLACEHOLDER();
  printf("diff = %g\n", dznrmmax_dist(m, Av, slots));
  for (unsigned int i=0; i<slots; i++)
    ASSERT(cabs(m[i]-Av[i])<1e-5);

TEST_DO("exit");
  he_free_sk(&sk);
  he_free_pk(&pk);
  for (unsigned int rot=0; rot<hectx.slots; rot++)
    he_free_pk(&rk[rot]);
  he_free_pt(&pt);
  he_free_ct(&ct_v);
  he_free_ct(&ct_Av);
  hectx_exit();
TEST_DONE();

TEST_END();
}

int main(int argc, char *argv[])
{
  if (argc != 2) {
    fprintf(stderr, "usage: %s <sk/pk>\n", argv[0]);
    exit(1);
  }
  MPI q=mpi_new(0);
  mpi_set_ui(q, 1);
  mpi_lshift(q, q, 220);
  unsigned int logn = 14;
  unsigned int slots = 8;
  uint64_t Delta = 1UL<<50;
  argc--;
  argv++;
  test_gemv(logn, q, slots, Delta, *argv);
  mpi_release(q);
  return 0;
}
