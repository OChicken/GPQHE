/*
 * Test file: homomorphic conjugate and rotation.
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

void blas_zrot(_Complex double y[], _Complex double x[], const unsigned int len, const int rot)
{
  assert(y!=x);
  for (int i=0; i<len; i++) {
    int idx = (i+rot)%len;
    if (idx<0)
      idx += len;
    assert(idx>=0);
    y[i] = x[idx];
  }
}

__attribute__((unused))
static void test_conj(const unsigned int logn,
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
  he_ct_t ct;
  he_alloc_ct(&ct);
  he_pk_t pk, ck;
  he_alloc_pk(&pk);
  he_alloc_pk(&ck);
  poly_mpi_t sk;
  he_alloc_sk(&sk);
  _Complex double m[slots], m0[slots];
  sample_z01vec(m0, slots);
  memcpy(m, m0, sizeof(m0));
  he_keypair(&pk, &sk);
TEST_DO("gen ck");
  he_genck(&ck, &sk);
TEST_DONE();

TEST_DO("conj(m)");
  for (unsigned int i=0; i<slots; i++)
    m0[i] = conj(m0[i]);
  he_ecd(&pt, m);

  if (!strcmp(key, "sk"))
    he_enc_sk(&ct, &pt, &sk);
  if (!strcmp(key, "pk"))
    he_enc_pk(&ct, &pt, &pk);
  he_show_ct_params(&ct, "ct (%s)", key);
  he_conj(&ct, &ck);
  he_show_ct_params(&ct, "conj(ct)");
  he_dec(&pt, &ct, &sk);
  he_dcd(m, &pt);
  TEST_PLACEHOLDER();
  printf("diff = %g\n", dznrmmax_dist(m, m0, slots));
  ASSERT(dznrmmax_dist(m, m0, slots)<1e-5);
  for (unsigned int i=0; i<slots; i++)
    ASSERT(cabs(m[i]-m0[i])<1e-5);
TEST_DONE();

TEST_DO("exit");
  he_free_sk(&sk);
  he_free_pk(&pk);
  he_free_pk(&ck);
  he_free_pt(&pt);
  he_free_ct(&ct);
  hectx_exit();
TEST_DONE();

TEST_END();
}

static void test_rot(const unsigned int logn,
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
  he_ct_t ct, ct0;
  he_alloc_ct(&ct);
  he_alloc_ct(&ct0);
  he_pk_t pk, rk[hectx.slots];
  he_alloc_pk(&pk);
  for (unsigned int rot=0; rot<hectx.slots; rot++)
    he_alloc_pk(&rk[rot]);
  poly_mpi_t sk;
  he_alloc_sk(&sk);
  _Complex double m[slots], m0[slots], mr[slots];
  sample_z01vec(m0, slots);
  he_keypair(&pk, &sk);
TEST_DO("gen rk");
  he_genrk(rk, &sk);
TEST_DONE();

TEST_DO("rot(m)");
  he_ecd(&pt, m0);
  if (!strcmp(key, "sk"))
    he_enc_sk(&ct0, &pt, &sk);
  if (!strcmp(key, "pk"))
    he_enc_pk(&ct0, &pt, &pk);

  for (unsigned int rot=0; rot<hectx.slots; rot++) {
    blas_zrot(mr, m0, slots, rot);
    he_copy_ct(&ct, &ct0);
    he_rot(&ct, rot, rk);
    he_dec(&pt, &ct, &sk);
    he_dcd(m, &pt);
    TEST_PLACEHOLDER();
    printf("diff[%d] = %g\n", rot, dznrmmax_dist(m, mr, slots));
    ASSERT(dznrmmax_dist(m, mr, slots)<1e-5);
  }
TEST_DONE();

TEST_DO("exit");
  he_free_sk(&sk);
  he_free_pk(&pk);
  for (unsigned int rot=0; rot<hectx.slots; rot++)
    he_free_pk(&rk[rot]);
  he_free_pt(&pt);
  he_free_ct(&ct);
  he_free_ct(&ct0);
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
  unsigned int slots = 64;
  uint64_t Delta = 1UL<<50;
  argc--;
  argv++;
  test_conj(logn, q, slots, Delta, *argv);
  test_rot (logn, q, slots, Delta, *argv);
  mpi_release(q);
  return 0;
}
