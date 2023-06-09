/*
 * Test file: poly mult check with self designed parameters.
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

void poly_print(poly_mpi_t *r)
{
  for (unsigned int k=0; k<polyctx.n; k++) {
    unsigned int i = polyctx.n-k-1;
    int64_t coeffs = mpi_to_s64(r->coeffs[i]);
    if (k) {
      if (coeffs>=0)
        printf(" + ");
      else
        printf(" - ");
    }
    printf("%li", __builtin_labs(coeffs));
    if ((i!=1)&&(i!=0))
      printf("*x^%u", i);
    else if (i==1)
      printf("*x");
  }
  printf("\n");
}

static void test_mult()
{
  poly_mpi_t a;
  poly_mpi_alloc(&a);
  poly_mpi_t b;
  poly_mpi_alloc(&b);
  poly_mpi_t r;
  poly_mpi_alloc(&r);

TEST_DO("a.coeffs[i] = i+2, b.coeffs[i] = i+3");
  for (unsigned int i=0; i<polyctx.n; i++) {
    mpi_set_ui(a.coeffs[i], i+2);
    mpi_set_ui(b.coeffs[i], i+3);
  }
  poly_mul(&r, &a, &b, polyctx.dimub, polyctx.q);
  poly_print(&r);
TEST_DONE();

TEST_DO("a.coeffs[i] = ctx.P[0]-i-1, b.coeffs[i] = ctx.P[1]-i-1");
  for (unsigned int i=0; i<polyctx.n; i++) {
    struct rns_ctx *rns=polyctx.rns;
    mpi_set_ui(a.coeffs[i], rns->p-i-1);
    rns=rns->next;
    mpi_set_ui(b.coeffs[i], rns->p-i-1);
  }
  poly_mul(&r, &a, &b, polyctx.dimub, polyctx.q);
  poly_print(&r);
TEST_DONE();

  poly_mpi_free(&a);
  poly_mpi_free(&b);
  poly_mpi_free(&r);
}

int main()
{
TEST_BEGIN();

  unsigned int logq = 61;
  MPI q = mpi_new(0);
  mpi_set_ui(q, 1);
  mpi_lshift(q, q, logq);
TEST_DO("init");
  polyctx_init(7,q);
TEST_DONE();

TEST_DO("check crt parameters");
  TEST_PLACEHOLDER();
  printf("logn = %u. n = %u\n", polyctx.logn, polyctx.n);
  TEST_PLACEHOLDER();
  printf("ctx.p = ");
  struct rns_ctx *rns = polyctx.rns;
  for (unsigned int d=0; d<polyctx.dimub; d++, rns=rns->next)
    printf("%lu ", rns->p);
  printf("\n");

  for (struct rns_ctx *rns=polyctx.rns; rns; rns=rns->next){
    TEST_PLACEHOLDER();
    printf("ctx.phat_invmp[%d] = ", rns->dim-1);
    for (unsigned int d=0; d<rns->dim; d++)
      printf("%lu ", rns->phat_invmp[d]);
    printf("\n");
  }
TEST_DONE();

  test_mult();

TEST_DO("exit");
  polyctx_exit();
TEST_DONE();

  gcry_mpi_release(q);
TEST_END();
  return 0;
}
