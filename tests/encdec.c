/*
 * Test file: encrypt and decrypt.
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

static void test_encdec(const unsigned int logn,
                        const MPI q,
                        const unsigned int slots,
                        const uint64_t Delta,
                        const char *key)
{
TEST_BEGIN();

TEST_DO("init");
  hectx_init(logn, q, slots, Delta);
TEST_DONE();
  he_pt_t pt; he_alloc_pt(&pt);
  he_pk_t pk; he_alloc_pk(&pk);
  he_ct_t ct; he_alloc_ct(&ct);
  poly_mpi_t sk;   he_alloc_sk(&sk);
  _Complex double m[slots], m0[slots];
  sample_z01vec(m, slots);
  memcpy(m0, m, sizeof(m));

TEST_DO("keypair");
  he_keypair(&pk, &sk);
TEST_DONE();

TEST_DO("encode");
  he_ecd(&pt, m);
TEST_DONE();
  he_show_pt_params(&pt, "pt before encrypt");

TEST_DO("encrypt");
  if (!strcmp(key, "sk"))
    he_enc_sk(&ct, &pt, &sk);
  if (!strcmp(key, "pk"))
    he_enc_pk(&ct, &pt, &pk);
TEST_DONE();
  he_show_ct_params(&ct, "ct encrypted by %s", key);

TEST_DO("decrypt");
  he_dec(&pt, &ct, &sk);
TEST_DONE();
  he_show_pt_params(&pt, "pt after decrypt");

TEST_DO("decode");
  he_dcd(m, &pt);
  double diff = dznrmmax_dist(m, m0, slots);
  TEST_PLACEHOLDER();
  printf("diff = %g\n", diff);
  for (unsigned int i=0; i<slots; i++)
    ASSERT(cabs(m[i]-m0[i])<1e-5);
TEST_DONE();

TEST_DO("exit");
  he_free_sk(&sk);
  he_free_pk(&pk);
  he_free_pt(&pt);
  he_free_ct(&ct);
  hectx_exit();
TEST_DONE();

TEST_END();
}

int main(int argc, const char *argv[])
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
  test_encdec(logn, q, slots, Delta, *argv);
  mpi_release(q);
  return 0;
}
