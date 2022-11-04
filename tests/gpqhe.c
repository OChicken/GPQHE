/*
 * Test file: homomorphic encryption, decryption, addition and multiplication.
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

/* canemb.c */
extern void canemb(_Complex double a[], const unsigned int slots);
extern void invcanemb(_Complex double a[], const unsigned int slots);

/* encode.c */
extern double canemb_norm(const _Complex double m[], const double Delta);
extern double canemb_norm_(const poly_mpi_t *a);

static inline double
blas_dznrmmax(const _Complex double z[], const unsigned int len)
{
  double norm = 0;
  for (unsigned int i=0; i<len; i++) {
    double tmp = cabs(z[i]);
    if (norm<tmp)
      norm = tmp;
  }
  return norm;
}

static inline double
blas_dnrmmax_dist(const double y[], const double x[], const unsigned int len)
{
  double dist = 0;
  for (unsigned int i=0; i<len; i++) {
    double tmp = fabs(y[i] - x[i]);
    if (dist<tmp)
      dist = tmp;
  }
  return dist;
}

static inline double
blas_dznrmmax_dist(const _Complex double y[], const _Complex double x[], const unsigned int len)
{
  double dist = 0;
  for (unsigned int i=0; i<len; i++) {
    double tmp = cabs(y[i] - x[i]);
    if (dist<tmp)
      dist = tmp;
  }
  return dist;
}

static inline void
blas_dzrot(_Complex double y[], const _Complex double x[], const unsigned int len, const int rot)
{
  assert(y!=x);
  for (unsigned int i=0; i<len; i++) {
    int idx = ((int)i+rot)%len;
    if (idx<0)
      idx += len;
    assert(idx>=0);
    y[i] = x[idx];
  }
}

static unsigned int slots;
static he_pk_t pk;
static poly_mpi_t sk;
static double diff;

#define CHECK_DIFF(m,m0)                       \
  TEST_PLACEHOLDER();                          \
  diff = blas_dznrmmax_dist((m), (m0), slots); \
  printf("diff = %g\n", diff);                 \
  ASSERT(diff<1e-5);

static void test_ecd()
{
TEST_BEGIN();
  he_pt_t pt, pt1, pt2;
  he_alloc_pt(&pt);
  he_alloc_pt(&pt1);
  he_alloc_pt(&pt2);
  _Complex double m[slots], m0[slots], m1[slots], m2[slots];

TEST_DO("canemb correctness");
  sample_z01vec(m0, slots);
  memcpy(m, m0, sizeof(m0));
  invcanemb(m,slots);
  canemb(m,slots);
  CHECK_DIFF(m, m0);
TEST_DONE();

TEST_DO("encode/decode");
  memcpy(m, m0, sizeof(m0));
  he_ecd(&pt, m);
  he_dcd(m, &pt);
  CHECK_DIFF(m, m0);
TEST_DONE();

TEST_DO("canemb_norm");
  sample_z01vec(m1, slots);
  sample_z01vec(m2, slots);
  for (unsigned int i=0; i<slots; i++) {
    m1[i] = (1-(double)i/slots)+(1-(double)i/slots)*I;
    m2[i] = (0.99-(double)i/slots)+(0.99-(double)i/slots)*I;
  }
  for (unsigned int i=0; i<slots; i++)
    m0[i] = m1[i]*m2[i];
  he_ecd(&pt1, m1);
  he_ecd(&pt2, m2);
  /* 按照论文的意思, Delta, 即nu, 应该大于等于canemb_norm */
  /* 实际上下列语句的输出结果取决于输入数据. 若输入数据小于1(可通过归一化手段), 则几乎肯定是大于的 */
  if (hectx.Delta>canemb_norm(m1, hectx.Delta))
    printf("Delta > canemb_norm_(m1, Delta)\n");
  else
    printf("Delta < canemb_norm_(m1, Delta)\n");
  if (hectx.Delta>canemb_norm_(&pt1.m))
    printf("Delta > canemb_norm(pt1)\n");
  else
    printf("Delta < canemb_norm(pt1)\n");
  if (hectx.Delta>canemb_norm(m2, hectx.Delta))
    printf("Delta > canemb_norm_(m2, Delta)\n");
  else
    printf("Delta < canemb_norm_(m2, Delta)\n");
  if (hectx.Delta>canemb_norm_(&pt2.m))
    printf("Delta > canemb_norm(pt2)\n");
  else
    printf("Delta < canemb_norm(pt2)\n");
  printf("pt1.nu                    = %f\n", pt1.nu);
  printf("canemb_norm(m1, Delta)    = %f\n", canemb_norm(m1, hectx.Delta));
  printf("canemb_norm(pt1)          = %f\n", canemb_norm_(&pt1.m));
  printf("pt2.nu                    = %f\n", pt2.nu);
  printf("canemb_norm(m2, Delta)    = %f\n", canemb_norm(m2, hectx.Delta));
  printf("canemb_norm(pt2)          = %f\n", canemb_norm_(&pt2.m));
  unsigned int dim = (mpi_get_nbits(hectx.q[hectx.L])*2)/GPQHE_LOGP+1;
  poly_mul(&pt.m, &pt1.m, &pt2.m, dim, hectx.q[hectx.L]);
  pt.nu = pt1.nu*pt2.nu;
  printf("pt.nu                     = %f\n", pt.nu);
  printf("canemb_norm_(m0, Delta^2) = %f\n", canemb_norm(m0, (double)hectx.Delta*hectx.Delta));
  printf("canemb_norm(pt)           = %f\n", canemb_norm_(&pt.m));
  he_dcd(m, &pt);
  CHECK_DIFF(m,m0);
TEST_DONE();

  he_free_pt(&pt);
  he_free_pt(&pt1);
  he_free_pt(&pt2);
TEST_END();
}

static void test_enc(const char *key)
{
TEST_BEGIN();
  he_pt_t pt; he_alloc_pt(&pt);
  he_pk_t pk; he_alloc_pk(&pk);
  he_ct_t ct; he_alloc_ct(&ct);
  poly_mpi_t sk;   he_alloc_sk(&sk);
  _Complex double m[slots], m0[slots];
  sample_z01vec(m, slots);
  memcpy(m0, m, sizeof(m));

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
TEST_DONE();
  CHECK_DIFF(m,m0);

TEST_DO("moddown");
  he_moddown(&ct);
  he_dec(&pt, &ct, &sk);
  he_dcd(m, &pt);
  CHECK_DIFF(m,m0);
TEST_DONE();

TEST_DO("exit");
  he_free_pt(&pt);
  he_free_ct(&ct);
TEST_DONE();

TEST_END();
}

static void test_add(const char *key)
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
    he_enc_sk(&ct1, &pt1, &sk);
    he_enc_sk(&ct2, &pt2, &sk);
  }
  if (!strcmp(key, "pk")) {
    he_enc_pk(&ct1, &pt1, &pk);
    he_enc_pk(&ct2, &pt2, &pk);
  }
  he_add(&ct, &ct1, &ct2);
  he_show_ct_params(&ct1, "ct1 (%s)", key);
  he_show_ct_params(&ct2, "ct2 (%s)", key);
  he_show_ct_params(&ct, "ct=ct1+ct2");
  he_dec(&pt, &ct, &sk);
  he_dcd(m, &pt);
  CHECK_DIFF(m,m0);
TEST_DONE();

TEST_DO("enc(m1)+m2");
  for (unsigned int i=0; i<slots; i++)
    m0[i] = m1[i]+m2[i];
  he_ecd(&pt1, m1);
  he_ecd(&pt2, m2);
  if (!strcmp(key, "sk"))
    he_enc_sk(&ct1, &pt1, &sk);
  if (!strcmp(key, "pk"))
    he_enc_pk(&ct1, &pt1, &pk);
  he_addpt(&ct, &ct1, &pt2);
  he_show_ct_params(&ct1, "ct1 (%s)", key);
  he_show_pt_params(&pt2, "pt2");
  he_show_ct_params(&ct, "ct=ct1+pt2");
  he_dec(&pt, &ct, &sk);
  he_dcd(m, &pt);
  CHECK_DIFF(m,m0);
TEST_DONE();

TEST_DO("enc(m1)-enc(m2)");
  for (unsigned int i=0; i<slots; i++)
    m0[i] = m1[i]-m2[i];
  he_ecd(&pt1, m1);
  he_ecd(&pt2, m2);
  if (!strcmp(key, "sk")) {
    he_enc_sk(&ct1, &pt1, &sk);
    he_enc_sk(&ct2, &pt2, &sk);
  }
  if (!strcmp(key, "pk")) {
    he_enc_pk(&ct1, &pt1, &pk);
    he_enc_pk(&ct2, &pt2, &pk);
  }
  he_sub(&ct, &ct1, &ct2);
  he_show_ct_params(&ct1, "ct1 (%s)", key);
  he_show_ct_params(&ct2, "ct2 (%s)", key);
  he_show_ct_params(&ct, "ct=ct1-ct2");
  he_dec(&pt, &ct, &sk);
  he_dcd(m, &pt);
  CHECK_DIFF(m,m0);
TEST_DONE();

TEST_DO("enc(m1)-m2");
  for (unsigned int i=0; i<slots; i++)
    m0[i] = m1[i]-m2[i];
  he_ecd(&pt1, m1);
  he_ecd(&pt2, m2);
  if (!strcmp(key, "sk"))
    he_enc_sk(&ct1, &pt1, &sk);
  if (!strcmp(key, "pk"))
    he_enc_pk(&ct1, &pt1, &pk);
  he_subpt(&ct, &ct1, &pt2);
  he_show_ct_params(&ct1, "ct1 (%s)", key);
  he_show_pt_params(&pt2, "pt2");
  he_show_ct_params(&ct, "ct=ct1-pt2");
  he_dec(&pt, &ct, &sk);
  he_dcd(m, &pt);
  CHECK_DIFF(m,m0);
TEST_DONE();

TEST_DO("m1-enc(m2)");
  for (unsigned int i=0; i<slots; i++)
    m0[i] = m1[i]-m2[i];
  he_ecd(&pt1, m1);
  he_ecd(&pt2, m2);
  if (!strcmp(key, "sk"))
    he_enc_sk(&ct2, &pt2, &sk);
  if (!strcmp(key, "pk"))
    he_enc_pk(&ct2, &pt2, &pk);
  he_subpt(&ct, &ct2, &pt1);
  he_neg(&ct);
  he_show_pt_params(&pt1, "pt1");
  he_show_ct_params(&ct2, "ct2 (%s)", key);
  he_show_ct_params(&ct, "ct=pt1-ct2");
  he_dec(&pt, &ct, &sk);
  he_dcd(m, &pt);
  CHECK_DIFF(m,m0);
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

static void test_mul(const char *key)
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
  he_evk_t rlk;
  he_alloc_evk(&rlk);
  _Complex double m[slots], m0[slots], m1[slots], m2[slots];
  sample_z01vec(m1, slots);
  sample_z01vec(m2, slots);

  for (unsigned int i=0; i<slots; i++)
    m0[i] = m1[i]*m2[i];
  he_ecd(&pt1, m1);
  he_ecd(&pt2, m2);
  if (!strcmp(key, "sk")) {
    he_enc_sk(&ct1, &pt1, &sk);
    he_enc_sk(&ct2, &pt2, &sk);
  }
  if (!strcmp(key, "pk")) {
    he_enc_pk(&ct1, &pt1, &pk);
    he_enc_pk(&ct2, &pt2, &pk);
  }
  he_show_ct_params(&ct1, "ct1 (%s)", key);
  he_show_ct_params(&ct2, "ct2 (%s)", key);

TEST_DO("gen rlk");
  he_genrlk(&rlk, &sk);
TEST_DONE();

TEST_DO("enc(m1)*enc(m2) without rescale");
  he_mul(&ct, &ct1, &ct2, &rlk);
  he_show_ct_params(&ct, "ct=ct1*ct2");
  he_dec(&pt, &ct, &sk);
  he_dcd(m, &pt);
  CHECK_DIFF(m,m0);
TEST_DONE();

TEST_DO("enc(m1)*enc(m2)");
  he_mul(&ct, &ct1, &ct2, &rlk);
  he_show_ct_params(&ct, "ct=ct1*ct2");
  he_rs(&ct);
  he_show_ct_params(&ct, "ct=ct1*ct2");
  he_dec(&pt, &ct, &sk);
  he_dcd(m, &pt);
  CHECK_DIFF(m,m0);
TEST_DONE();

TEST_DO("enc(m1)*enc(m2) moddown correctness");
  he_mul(&ct, &ct1, &ct2, &rlk);
  he_rs(&ct);
  he_moddown(&ct);
  he_show_ct_params(&ct, "ct=ct1*ct2");
  he_dec(&pt, &ct, &sk);
  he_dcd(m, &pt);
  CHECK_DIFF(m,m0);
TEST_DONE();

TEST_DO("enc(m1)*enc(m2) moddown correctness with ct1 ct2 are moddown");
  he_moddown(&ct1);
  he_moddown(&ct2);
  he_show_ct_params(&ct1, "ct1 (%s)", key);
  he_show_ct_params(&ct2, "ct2 (%s)", key);
  he_mul(&ct, &ct1, &ct2, &rlk);
  he_rs(&ct);
  he_show_ct_params(&ct, "ct=ct1*ct2");
  he_dec(&pt, &ct, &sk);
  he_dcd(m, &pt);
  CHECK_DIFF(m,m0);
TEST_DONE();

  if (!strcmp(key, "sk"))
    he_enc_sk(&ct1, &pt1, &sk);
  if (!strcmp(key, "pk"))
    he_enc_pk(&ct1, &pt1, &pk);

TEST_DO("enc(m1)*m2");
  he_mulpt(&ct, &ct1, &pt2);
  he_rs(&ct);
  he_show_ct_params(&ct1, "ct1 (%s)", key);
  he_show_pt_params(&pt2, "pt2");
  he_show_ct_params(&ct, "ct=ct1*pt2");
  he_dec(&pt, &ct, &sk);
  he_dcd(m, &pt);
  CHECK_DIFF(m,m0);
TEST_DONE();

TEST_DO("enc(m1)*m2 with ct1 is moddown");
  he_moddown(&ct1);
  he_mulpt(&ct, &ct1, &pt2);
  he_rs(&ct);
  he_show_ct_params(&ct1, "ct1 (%s)", key);
  he_show_pt_params(&pt2, "pt2");
  he_show_ct_params(&ct, "ct=ct1*ct2");
  he_dec(&pt, &ct, &sk);
  he_dcd(m, &pt);
  CHECK_DIFF(m,m0);
TEST_DONE();

TEST_DO("exit");
  he_free_evk(&rlk);
  he_free_pt(&pt);
  he_free_pt(&pt1);
  he_free_pt(&pt2);
  he_free_ct(&ct);
  he_free_ct(&ct1);
  he_free_ct(&ct2);
TEST_DONE();

TEST_END();
}

static void test_conj(const char *key)
{
TEST_BEGIN();
  he_pt_t pt;
  he_alloc_pt(&pt);
  he_ct_t ct, ct0;
  he_alloc_ct(&ct);
  he_alloc_ct(&ct0);
  he_evk_t ck;
  he_alloc_evk(&ck);
  _Complex double m[slots], m0[slots];
  sample_z01vec(m0, slots);
  he_ecd(&pt, m0);
  if (!strcmp(key, "sk"))
    he_enc_sk(&ct0, &pt, &sk);
  if (!strcmp(key, "pk"))
    he_enc_pk(&ct0, &pt, &pk);
  he_show_ct_params(&ct0, "ct0 (%s)", key);
  memcpy(m, m0, sizeof(m0));
  for (unsigned int i=0; i<slots; i++)
    m0[i] = conj(m0[i]);

TEST_DO("gen ck");
  he_genck(&ck, &sk);
TEST_DONE();

  he_copy_ct(&ct, &ct0);
TEST_DO("conj(m)");
  he_conj(&ct, &ck);
  he_show_ct_params(&ct, "conj(ct)");
  he_dec(&pt, &ct, &sk);
  he_dcd(m, &pt);
  CHECK_DIFF(m,m0);
TEST_DONE();

  he_copy_ct(&ct, &ct0);
  he_moddown(&ct);
  he_moddown(&ct);
TEST_DO("conj(m) with ckhat with ct has been moddown");
  he_conj(&ct, &ck);
  he_show_ct_params(&ct, "conj(ct)");
  he_dec(&pt, &ct, &sk);
  he_dcd(m, &pt);
  CHECK_DIFF(m,m0);
TEST_DONE();

TEST_DO("exit");
  he_free_evk(&ck);
  he_free_pt(&pt);
  he_free_ct(&ct);
  he_free_ct(&ct0);
TEST_DONE();

TEST_END();
}

static void test_rot(const char *key)
{
TEST_BEGIN();
  he_pt_t pt;
  he_alloc_pt(&pt);
  he_ct_t ct, ct0;
  he_alloc_ct(&ct);
  he_alloc_ct(&ct0);
  he_evk_t rk[hectx.slots];
  for (unsigned int rot=0; rot<hectx.slots; rot++)
    he_alloc_evk(&rk[rot]);
  _Complex double m[slots], m0[slots], mr[slots];
  sample_z01vec(m0, slots);
  he_ecd(&pt, m0);
  if (!strcmp(key, "sk"))
    he_enc_sk(&ct0, &pt, &sk);
  if (!strcmp(key, "pk"))
    he_enc_pk(&ct0, &pt, &pk);
  he_show_ct_params(&ct0, "ct0 (%s)", key);

TEST_DO("gen rk");
  he_genrk(rk, &sk);
TEST_DONE();

TEST_DO("rot(m)");
  for (unsigned int rot=0; rot<hectx.slots; rot++) {
    blas_dzrot(mr, m0, slots, rot);
    he_copy_ct(&ct, &ct0);
    he_rot(&ct, rot, rk);
    he_dec(&pt, &ct, &sk);
    he_dcd(m, &pt);
    TEST_PLACEHOLDER();
    diff = blas_dznrmmax_dist(m, mr, slots);
    printf("diff[%d] = %g\n", rot, diff);
    ASSERT(diff<1e-5);
  }
TEST_DONE();

TEST_DO("exit");
  for (unsigned int rot=0; rot<hectx.slots; rot++)
    he_free_evk(&rk[rot]);
  he_free_pt(&pt);
  he_free_ct(&ct);
  he_free_ct(&ct0);
TEST_DONE();

TEST_END();
}

static void test_gemv(const char *key)
{
TEST_BEGIN();
  he_pt_t pt;
  he_alloc_pt(&pt);
  he_ct_t ct_Av, ct_v;
  he_alloc_ct(&ct_Av);
  he_alloc_ct(&ct_v);
  he_evk_t rk[hectx.slots];
  for (unsigned int rot=0; rot<hectx.slots; rot++)
    he_alloc_evk(&rk[rot]);
  _Complex double m[slots], Av[slots], A[slots*slots], v[slots];
  sample_z01vec(v, slots);
  sample_z01vec(A, slots*slots);
  memset(Av, 0, sizeof(Av));
  for (unsigned int i=0; i<slots; i++)
    for (unsigned int j=0; j<slots; j++)
      Av[i] += A[i*slots+j]*v[j];
  he_ecd(&pt, v);
  if (!strcmp(key, "sk"))
    he_enc_sk(&ct_v, &pt, &sk);
  if (!strcmp(key, "pk"))
    he_enc_pk(&ct_v, &pt, &pk);

TEST_DO("gen rk");
  he_genrk(rk, &sk);
TEST_DONE();

TEST_DO("gemv");
  he_gemv(&ct_Av, A, &ct_v, rk);
  he_dec(&pt, &ct_Av, &sk);
  he_dcd(m, &pt);
  CHECK_DIFF(m, Av);
TEST_DONE();

TEST_DO("exit");
  for (unsigned int rot=0; rot<hectx.slots; rot++)
    he_free_evk(&rk[rot]);
  he_free_pt(&pt);
  he_free_ct(&ct_v);
  he_free_ct(&ct_Av);
TEST_DONE();

TEST_END();
}

static void test_inv(const char *key, const unsigned int iter)
{
TEST_BEGIN();
  he_pt_t pt;
  he_alloc_pt(&pt);
  he_ct_t ct_inv, ct;
  he_alloc_ct(&ct_inv);
  he_alloc_ct(&ct);
  he_evk_t rlk;
  he_alloc_evk(&rlk);
  _Complex double m[slots], mraw[slots], m0[slots], minv[slots], an[slots], bn[slots];
  sample_z01vec(m0, slots);
  sample_z01vec(m0, slots);
  //sample_z01vec(m0, slots);
  for (unsigned int i=0; i<slots; i++) {
    m0[i] = creal(m0[i]);
    minv[i] = 1/m0[i];
  }

TEST_DO("gen rlk");
  he_genrlk(&rlk, &sk);
TEST_DONE();

  memcpy(m, m0, sizeof(m0));
TEST_DO("raw");
  for (unsigned int i=0; i<slots; i++) {
    an[i] = 2-m[i];
    bn[i] = 1-m[i];
    for (unsigned int k=0; k<iter; k++) {
      bn[i] = bn[i]*bn[i];
      an[i] = an[i]*(bn[i]+1);
    }
    m[i] = an[i];
  }
  CHECK_DIFF(m, minv);
TEST_DONE();

  memcpy(m, m0, sizeof(m0));
  he_ecd(&pt, m);
  if (!strcmp(key, "sk"))
    he_enc_sk(&ct, &pt, &sk);
  if (!strcmp(key, "pk"))
    he_enc_pk(&ct, &pt, &pk);

TEST_DO("inv");
  he_inv(&ct_inv, &ct, &rlk, iter);
  he_show_ct_params(&ct_inv, "ct_inv");
  he_dec(&pt, &ct_inv, &sk);
  he_dcd(m, &pt);
  CHECK_DIFF(m, minv);
TEST_DONE();

TEST_DO("exit");
  he_free_pt(&pt);
  he_free_ct(&ct);
  he_free_ct(&ct_inv);
  he_free_evk(&rlk);
TEST_DONE();

TEST_END();
}

void set(unsigned int *logn, unsigned int *logq, unsigned int *slots,
  unsigned int *logDelta, unsigned int *iter, int argc, char *argv[])
{
  if ((argc==1) || (argc==2&&strcmp(argv[1], "ecd"))) {
    fprintf(stderr, "usage: %s [ecd/enc/add/mul/conj/rot/gemv/inv] [sk/pk] "
      "--logn=num --logq=num --logp=num --iter=num\n", argv[0]);
    exit(1);
  }
  unsigned int idx=3; /* start index of parameter chosen */
  while (argv[idx]) {
    if (!strncmp(argv[idx], "--logn=", 7))
      *logn=atoi(argv[idx]+7), printf("argv[%u]=%s %u\n", idx, argv[idx], __LINE__), idx++;
    else if (!strncmp(argv[idx], "--logq=", 7))
      *logq=atoi(argv[idx]+7), idx++;
    else if (!strncmp(argv[idx], "--slots=", 8))
      *slots=atoi(argv[idx]+8), idx++;
    else if (!strncmp(argv[idx], "--logDelta=", 11))
      *logDelta=atoi(argv[idx]+11), idx++;
    else if (!strncmp(argv[idx], "--iter=", 7))
      *iter=atoi(argv[idx]+7), idx++;
  }
}

int main(int argc, char *argv[])
{
  MPI q = mpi_set_ui(NULL, 1);
  unsigned int logn = 14;
  unsigned int logq = 220;
  unsigned int logDelta = 50;
  unsigned int iter;
  slots = 64;
  if (!strcmp(argv[1], "inv")) {
    logq  = 438;
    slots = 4;
    iter  = 6;
  }
  set(&logn, &logq, &slots, &logDelta, &iter, argc, argv);
  mpi_lshift(q, q, logq);
  uint64_t Delta = 1UL<<logDelta;
  hectx_init(logn, q, slots, Delta);
  he_show_ctx_params();
  /* init keys */
  he_alloc_pk(&pk);
  he_alloc_sk(&sk);
  /* main */
  if (!strcmp(argv[1], "ecd")){
    test_ecd();
    goto release;
  }
  he_keypair(&pk, &sk);
  if (!strcmp(argv[1], "enc"))
    test_enc (argv[2]);
  if (!strcmp(argv[1], "add"))
    test_add (argv[2]);
  if (!strcmp(argv[1], "mul"))
    test_mul (argv[2]);
  if (!strcmp(argv[1], "conj"))
    test_conj(argv[2]);
  if (!strcmp(argv[1], "rot"))
    test_rot (argv[2]);
  if (!strcmp(argv[1], "gemv"))
    test_gemv(argv[2]);
  if (!strcmp(argv[1], "inv"))
    test_inv(argv[2], iter);

release:
  /* release */
  mpi_release(q);
  he_free_sk(&sk);
  he_free_pk(&pk);
  hectx_exit();
  return 0;
}
