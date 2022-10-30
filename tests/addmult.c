/*
 * Test file: homomorphic addition and multiplication.
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

static void test_addsub(const unsigned int slots,
                        const char *key,
                        const he_pk_t *pk,
                        const poly_mpi_t *sk)
{
TEST_BEGIN();
  he_pt_t pt, pt1, pt2;
  he_alloc_pt(&pt);
  he_alloc_pt(&pt1);
  he_alloc_pt(&pt2);
  he_ct_t ct, ct1, ct2;
  he_alloc_ct(&ct);
  he_alloc_ct(&ct1);
  he_alloc_ct(&ct2);

  _Complex double m[slots], m0[slots], m1[slots], m2[slots];
  sample_z01vec(m1, slots);
  sample_z01vec(m2, slots);

TEST_DO("enc(m1)+enc(m2)");
  for (unsigned int i=0; i<slots; i++)
    m0[i] = m1[i]+m2[i];
  he_ecd(&pt1, m1);
  he_ecd(&pt2, m2);
  if (!strcmp(key, "sk")) {
    he_enc_sk(&ct1, &pt1, sk);
    he_enc_sk(&ct2, &pt2, sk);
  }
  if (!strcmp(key, "pk")) {
    he_enc_pk(&ct1, &pt1, pk);
    he_enc_pk(&ct2, &pt2, pk);
  }
  he_add(&ct, &ct1, &ct2);
  he_show_ct_params(&ct1, "ct1 (%s)", key);
  he_show_ct_params(&ct2, "ct2 (%s)", key);
  he_show_ct_params(&ct, "ct=ct1+ct2");
  he_dec(&pt, &ct, sk);
  he_dcd(m, &pt);
  TEST_PLACEHOLDER();
  printf("diff = %g\n", dznrmmax_dist(m, m0, slots));
  for (unsigned int i=0; i<slots; i++)
    ASSERT(cabs(m[i]-m0[i])<1e-5);
TEST_DONE();

TEST_DO("enc(m1)+m2");
  for (unsigned int i=0; i<slots; i++)
    m0[i] = m1[i]+m2[i];
  he_ecd(&pt1, m1);
  he_ecd(&pt2, m2);
  if (!strcmp(key, "sk"))
    he_enc_sk(&ct1, &pt1, sk);
  if (!strcmp(key, "pk"))
    he_enc_pk(&ct1, &pt1, pk);
  he_addpt(&ct, &ct1, &pt2);
  he_show_ct_params(&ct1, "ct1 (%s)", key);
  he_show_pt_params(&pt2, "pt2");
  he_show_ct_params(&ct, "ct=ct1+pt2");
  he_dec(&pt, &ct, sk);
  he_dcd(m, &pt);
  TEST_PLACEHOLDER();
  printf("diff = %g\n", dznrmmax_dist(m, m0, slots));
  for (unsigned int i=0; i<slots; i++)
    ASSERT(cabs(m[i]-m0[i])<1e-5);
TEST_DONE();

TEST_DO("enc(m1)-enc(m2)");
  for (unsigned int i=0; i<slots; i++)
    m0[i] = m1[i]-m2[i];
  he_ecd(&pt1, m1);
  he_ecd(&pt2, m2);
  if (!strcmp(key, "sk")) {
    he_enc_sk(&ct1, &pt1, sk);
    he_enc_sk(&ct2, &pt2, sk);
  }
  if (!strcmp(key, "pk")) {
    he_enc_pk(&ct1, &pt1, pk);
    he_enc_pk(&ct2, &pt2, pk);
  }
  he_sub(&ct, &ct1, &ct2);
  he_show_ct_params(&ct1, "ct1 (%s)", key);
  he_show_ct_params(&ct2, "ct2 (%s)", key);
  he_show_ct_params(&ct, "ct=ct1-ct2");
  he_dec(&pt, &ct, sk);
  he_dcd(m, &pt);
  TEST_PLACEHOLDER();
  printf("diff = %g\n", dznrmmax_dist(m, m0, slots));
  for (unsigned int i=0; i<slots; i++)
    ASSERT(cabs(m[i]-m0[i])<1e-5);
TEST_DONE();

TEST_DO("enc(m1)-m2");
  for (unsigned int i=0; i<slots; i++)
    m0[i] = m1[i]-m2[i];
  he_ecd(&pt1, m1);
  he_ecd(&pt2, m2);
  if (!strcmp(key, "sk"))
    he_enc_sk(&ct1, &pt1, sk);
  if (!strcmp(key, "pk"))
    he_enc_pk(&ct1, &pt1, pk);
  he_subpt(&ct, &ct1, &pt2);
  he_show_ct_params(&ct1, "ct1 (%s)", key);
  he_show_pt_params(&pt2, "pt2");
  he_show_ct_params(&ct, "ct=ct1-pt2");
  he_dec(&pt, &ct, sk);
  he_dcd(m, &pt);
  TEST_PLACEHOLDER();
  printf("diff = %g\n", dznrmmax_dist(m, m0, slots));
  for (unsigned int i=0; i<slots; i++)
    ASSERT(cabs(m[i]-m0[i])<1e-5);
TEST_DONE();

TEST_DO("m1-enc(m2)");
  for (unsigned int i=0; i<slots; i++)
    m0[i] = m1[i]-m2[i];
  he_ecd(&pt1, m1);
  he_ecd(&pt2, m2);
  if (!strcmp(key, "sk"))
    he_enc_sk(&ct2, &pt2, sk);
  if (!strcmp(key, "pk"))
    he_enc_pk(&ct2, &pt2, pk);
  he_subpt(&ct, &ct2, &pt1);
  he_neg(&ct);
  he_show_pt_params(&pt1, "pt1");
  he_show_ct_params(&ct2, "ct2 (%s)", key);
  he_show_ct_params(&ct, "ct=pt1-ct2");
  he_dec(&pt, &ct, sk);
  he_dcd(m, &pt);
  TEST_PLACEHOLDER();
  printf("diff = %g\n", dznrmmax_dist(m, m0, slots));
  for (unsigned int i=0; i<slots; i++)
    ASSERT(cabs(m[i]-m0[i])<1e-5);
TEST_DONE();

TEST_DO("exit");
  he_free_pt(&pt);
  he_free_pt(&pt1);
  he_free_pt(&pt2);
  he_free_ct(&ct);
  he_free_ct(&ct1);
  he_free_ct(&ct2);
TEST_DONE();

TEST_END();
}

static void test_mult(const unsigned int slots,
                        const char *key,
                        const he_pk_t *pk,
                        const poly_mpi_t *sk)
{
TEST_BEGIN();

  struct he_pt pt, pt1, pt2;
  he_alloc_pt(&pt);
  he_alloc_pt(&pt1);
  he_alloc_pt(&pt2);
  struct he_ct ct, ct1, ct2;
  he_alloc_ct(&ct);
  he_alloc_ct(&ct1);
  he_alloc_ct(&ct2);
  struct he_pk rlk;
  he_alloc_pk(&rlk);
  _Complex double m[slots], m0[slots], m1[slots], m2[slots];
  sample_z01vec(m1, slots);
  sample_z01vec(m2, slots);
  he_genrlk(&rlk, sk);

TEST_DO("enc(m1)*enc(m2)");
  for (unsigned int i=0; i<slots; i++)
    m0[i] = m1[i]*m2[i];
  he_ecd(&pt1, m1);
  he_ecd(&pt2, m2);
  if (!strcmp(key, "sk")) {
    he_enc_sk(&ct1, &pt1, sk);
    he_enc_sk(&ct2, &pt2, sk);
  }
  if (!strcmp(key, "pk")) {
    he_enc_pk(&ct1, &pt1, pk);
    he_enc_pk(&ct2, &pt2, pk);
  }
  he_show_ct_params(&ct1, "ct1 (%s)", key);
  he_show_ct_params(&ct2, "ct2 (%s)", key);
  he_mul(&ct, &ct1, &ct2, &rlk);
  he_show_ct_params(&ct, "ct=ct1*ct2");
  he_rs(&ct);
  he_show_ct_params(&ct, "ct=ct1*ct2");
  he_dec(&pt, &ct, sk);
  he_dcd(m, &pt);
  TEST_PLACEHOLDER();
  printf("diff = %g\n", dznrmmax_dist(m, m0, slots));
  for (unsigned int i=0; i<slots; i++)
    ASSERT(cabs(m[i]-m0[i])<1e-5);
TEST_DONE();

TEST_DO("enc(m1)*m2");
  for (unsigned int i=0; i<slots; i++)
    m0[i] = m1[i]*m2[i];
  he_ecd(&pt1, m1);
  he_ecd(&pt2, m2);
  if (!strcmp(key, "sk"))
    he_enc_sk(&ct1, &pt1, sk);
  if (!strcmp(key, "pk"))
    he_enc_pk(&ct1, &pt1, pk);
  he_mulpt(&ct, &ct1, &pt2);
  he_rs(&ct);
  he_show_ct_params(&ct1, "ct1 (%s)", key);
  he_show_pt_params(&pt2, "pt2");
  he_show_ct_params(&ct, "ct=ct1*pt2");
  he_dec(&pt, &ct, sk);
  he_dcd(m, &pt);
  TEST_PLACEHOLDER();
  printf("diff = %g\n", dznrmmax_dist(m, m0, slots));
  for (unsigned int i=0; i<slots; i++)
    ASSERT(cabs(m[i]-m0[i])<1e-5);
TEST_DONE();

TEST_DO("exit");

  he_free_pk(&rlk);
  he_free_pt(&pt);
  he_free_pt(&pt1);
  he_free_pt(&pt2);
  he_free_ct(&ct);
  he_free_ct(&ct1);
  he_free_ct(&ct2);
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
  hectx_init(logn, q, slots, Delta);
  he_show_ctx_params();
  /* init keys */
  he_pk_t pk;
  he_alloc_pk(&pk);
  poly_mpi_t sk;
  he_alloc_sk(&sk);
  he_keypair(&pk, &sk);
  /* main */
  argc--;
  argv++;
  test_addsub(slots, *argv, &pk, &sk);
  test_mult(slots, *argv, &pk, &sk);
  /* release */
  mpi_release(q);
  he_free_sk(&sk);
  he_free_pk(&pk);
  hectx_exit();
  return 0;
}
