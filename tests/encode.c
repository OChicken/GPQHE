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

/* poly.h */
extern struct poly_ctx polyctx;

/* canemb.c */
extern void canemb(_Complex double a[], const unsigned int slots);
extern void invcanemb(_Complex double a[], const unsigned int slots);

/* encode.c */
extern void encode(poly_mpi_t *r, const _Complex double m[], const double Delta);
extern void decode(_Complex double m[], poly_mpi_t *a, const double Delta, const MPI q);
extern double canemb_norm(const _Complex double m[], const double Delta);

static void test_encode()
{
TEST_BEGIN();

  MPI q=mpi_new(0);
  mpi_set_ui(q, 1);
  mpi_lshift(q, q, 220);
  unsigned int slots=8;
  uint64_t Delta = 1UL<<50;
TEST_DO("init");
  hectx_init(14, q, slots, Delta);
TEST_DONE();
  poly_mpi_t pt,pt1,pt2;
  poly_mpi_alloc(&pt);
  poly_mpi_alloc(&pt1);
  poly_mpi_alloc(&pt2);

  _Complex double m[slots],m0[slots],m1[slots],m2[slots];
  for (unsigned int i=0; i<slots; i++)
    m[i] = m0[i] = rand()+rand()*I;

TEST_DO("canemb");
  invcanemb(m,slots);
  canemb(m,slots);
  for (unsigned int i=0; i<slots; i++)
    assert(cabs(m[i]-m0[i])<1e-5);
TEST_DONE();

TEST_DO("encode/decode");
  TEST_PLACEHOLDER();
  printf("canemb norm = %.10g\n", canemb_norm(m, Delta));
  encode(&pt, m, Delta);
  decode(m, &pt, Delta, q);
  for (unsigned int i=0; i<slots; i++)
    assert(cabs(m[i]-m0[i])<1e-5);
TEST_DONE();

TEST_DO("conj");
  poly_mpi_t ptp;
  poly_mpi_alloc(&ptp);
  encode(&pt, m, Delta);
  poly_conj(&ptp, &pt);
  decode(m, &ptp, Delta, q);
  for (unsigned int i=0; i<slots; i++)
    assert(cabs(conj(m[i])-m0[i])<1e-5);
  poly_mpi_free(&ptp);
TEST_DONE();

TEST_DO("mult");
  for (unsigned int i=0; i<slots; i++) {
    m1[i] = (i+1)+(i+2)*I;
    m2[i] = (i+slots+1)+(i+slots+2)*I;
    m0[i] = m1[i]*m2[i];
  }
  encode(&pt1, m1, Delta);
  encode(&pt2, m2, Delta);
  unsigned int dim = mpi_get_nbits(polyctx.q)/GPQHE_LOGP+1;
  poly_mul(&pt, &pt1, &pt2, dim, q);
  decode(m, &pt, (double)Delta*(double)Delta, polyctx.q);
  for (unsigned int i=0; i<slots; i++) {
    ASSERT(cabs(m[i]-m0[i])<1e-5);
    TEST_PLACEHOLDER();
    printf("(%.0f,%.0f)*(%.0f,%.0f)=(%.0f,%.0f)≈(%.0f,%.0f)\n",
      creal(m1[i]), cimag(m1[i]),
      creal(m2[i]), cimag(m2[i]),
      creal(m0[i]), cimag(m0[i]),
      creal(m[i]), cimag(m[i]));
  }
TEST_DONE();

  poly_mpi_free(&pt);
  poly_mpi_free(&pt1);
  poly_mpi_free(&pt2);
  mpi_release(q);
  hectx_exit();
TEST_END();
}

int main()
{
  test_encode();
  return 0;
}
