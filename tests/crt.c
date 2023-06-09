/*
 * Test file: crt check with self designed parameters.
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

/* rns.c */
extern void rns_decompose(uint64_t ahat[], const MPI a[], const struct rns_ctx *rns);
extern void rns_reconstruct(MPI a[], const uint64_t ahat[], const unsigned int i, const struct rns_ctx *rns);

static void test_crt()
{
TEST_BEGIN();

  unsigned int logq = 10;
  MPI q = mpi_new(0);
  mpi_set_ui(q, 1);
  mpi_lshift(q, q, logq);
TEST_DO("init");
  polyctx_init(4,q);
TEST_DONE();

TEST_DO("check crt parameters");
  TEST_PLACEHOLDER();
  printf("ctx.p = ");
  for (struct rns_ctx *rns=polyctx.rns; rns; rns=rns->next)
    printf("%lu ", rns->p);
  printf("\n");

  TEST_PLACEHOLDER();
  printf("ctx.P = ");
  for (struct rns_ctx *rns=polyctx.rns; rns; rns=rns->next)
    printf("%lu ", mpi_to_u64(rns->P));
  printf("\n");

  for (struct rns_ctx *rns=polyctx.rns; rns; rns=rns->next) {
    TEST_PLACEHOLDER();
    printf("ctx.phat[%d] = ", rns->dim-1);
    for (unsigned int d=0; d<rns->dim; d++)
      printf("%lu ", mpi_to_u64(rns->phat[d]));
    printf("\n");
  }

  for (struct rns_ctx *rns=polyctx.rns; rns; rns=rns->next){
    TEST_PLACEHOLDER();
    printf("ctx.phat_invmp[%d] = ", rns->dim-1);
    for (unsigned int d=0; d<rns->dim; d++)
      printf("%lu ", rns->phat_invmp[d]);
    printf("\n");
  }
TEST_DONE();

  poly_mpi_t a;
  poly_mpi_alloc(&a);

TEST_DO("a.coeffs[i] = ctx.P[2]-i-1");
  struct rns_ctx *rns = polyctx.rns->next->next;
  for (unsigned int i=0; i<polyctx.n; i++)
    mpi_set_ui(a.coeffs[i], mpi_to_u64(rns->P)-i-1);
  poly_rns_t ahat;
  TEST_PLACEHOLDER();
  printf("a = ");
  for (unsigned int i=0; i<polyctx.n; i++)
    printf("%lu ", mpi_to_u64(a.coeffs[i]));
  printf("\n");
  MPI Q = mpi_new(0);
  rns = polyctx.rns->next->next->next->next->next;
  mpi_set(Q, rns->P);
TEST_DONE();

TEST_DO("rns decompose dim=dim_max=6, end index=5");
  unsigned int dim = polyctx.dimub;
  poly_rns_alloc(&ahat, dim);
  rns=polyctx.rns;
  for (unsigned int d=0; d<dim; d++) {
    rns_decompose(&ahat.coeffs[d*polyctx.n], a.coeffs, rns);
    rns=rns->next;
  }
  for (unsigned int d=0; d<dim; d++) {
    TEST_PLACEHOLDER();
    printf("ahat[%u] = ", d);
    for (unsigned int i=0; i<polyctx.n; i++)
      printf("%lu ", ahat.coeffs[d*polyctx.n+i]);
    printf("\n");
  }
  for (rns=polyctx.rns; rns->dim<dim; rns=rns->next);
  for (unsigned int i=0; i<polyctx.n; i++)
    rns_reconstruct(a.coeffs, ahat.coeffs, i, rns);
  TEST_PLACEHOLDER();
  printf("a = ");
  for (unsigned int i=0; i<polyctx.n; i++)
    printf("%lu ", mpi_to_u64(a.coeffs[i]));
  printf("\n");
  poly_rns_free(&ahat);
TEST_DONE();

TEST_DO("rns decompose dim=dim_max-1, end index=4");
  dim--;
  poly_rns_alloc(&ahat, dim);
  rns=polyctx.rns;
  for (unsigned int d=0; d<dim; d++, rns=rns->next)
    rns_decompose(&ahat.coeffs[d*polyctx.n], a.coeffs, rns);
  for (unsigned int d=0; d<dim; d++) {
    TEST_PLACEHOLDER();
    printf("ahat[%u] = ", d);
    for (unsigned int i=0; i<polyctx.n; i++)
      printf("%lu ", ahat.coeffs[d*polyctx.n+i]);
    printf("\n");
  }
  for (rns=polyctx.rns; rns->dim<dim; rns=rns->next);
  for (unsigned int i=0; i<polyctx.n; i++)
    rns_reconstruct(a.coeffs, ahat.coeffs, i, rns);
  TEST_PLACEHOLDER();
  printf("a = ");
  for (unsigned int i=0; i<polyctx.n; i++)
    printf("%lu ", mpi_to_u64(a.coeffs[i]));
  printf("\n");
  poly_rns_free(&ahat);
TEST_DONE();

TEST_DO("rns decompose dim=dim_max-2, end index=3");
  dim--;
  poly_rns_alloc(&ahat, dim);
  rns=polyctx.rns;
  for (unsigned int d=0; d<dim; d++, rns=rns->next)
    rns_decompose(&ahat.coeffs[d*polyctx.n], a.coeffs, rns);
  for (unsigned int d=0; d<dim; d++) {
    TEST_PLACEHOLDER();
    printf("ahat[%u] = ", d);
    for (unsigned int i=0; i<polyctx.n; i++)
      printf("%lu ", ahat.coeffs[d*polyctx.n+i]);
    printf("\n");
  }
  for (rns=polyctx.rns; rns->dim<dim; rns=rns->next);
  for (unsigned int i=0; i<polyctx.n; i++)
    rns_reconstruct(a.coeffs, ahat.coeffs, i, rns);
  TEST_PLACEHOLDER();
  printf("a = ");
  for (unsigned int i=0; i<polyctx.n; i++)
    printf("%lu ", mpi_to_u64(a.coeffs[i]));
  printf("\n");
  poly_rns_free(&ahat);
TEST_DONE();

TEST_DO("rns decompose dim=dim_max-3, end index=2");
  dim--;
  poly_rns_alloc(&ahat, dim);
  rns=polyctx.rns;
  for (unsigned int d=0; d<dim; d++, rns=rns->next)
    rns_decompose(&ahat.coeffs[d*polyctx.n], a.coeffs, rns);
  for (unsigned int d=0; d<dim; d++) {
    TEST_PLACEHOLDER();
    printf("ahat[%u] = ", d);
    for (unsigned int i=0; i<polyctx.n; i++)
      printf("%lu ", ahat.coeffs[d*polyctx.n+i]);
    printf("\n");
  }
  for (rns=polyctx.rns; rns->dim<dim; rns=rns->next);
  for (unsigned int i=0; i<polyctx.n; i++)
    rns_reconstruct(a.coeffs, ahat.coeffs, i, rns);
  TEST_PLACEHOLDER();
  printf("a = ");
  for (unsigned int i=0; i<polyctx.n; i++)
    printf("%lu ", mpi_to_u64(a.coeffs[i]));
  printf("\n");
  poly_rns_free(&ahat);
TEST_DONE();

TEST_DO("rns decompose dim=dim_max-4, end index=1");
  dim--;
  poly_rns_alloc(&ahat, dim);
  rns=polyctx.rns;
  for (unsigned int d=0; d<dim; d++, rns=rns->next)
    rns_decompose(&ahat.coeffs[d*polyctx.n], a.coeffs, rns);
  for (unsigned int d=0; d<dim; d++) {
    TEST_PLACEHOLDER();
    printf("ahat[%u] = ", d);
    for (unsigned int i=0; i<polyctx.n; i++)
      printf("%lu ", ahat.coeffs[d*polyctx.n+i]);
    printf("\n");
  }
  for (rns=polyctx.rns; rns->dim<dim; rns=rns->next);
  for (unsigned int i=0; i<polyctx.n; i++)
    rns_reconstruct(a.coeffs, ahat.coeffs, i, rns);
  TEST_PLACEHOLDER();
  printf("a = ");
  for (unsigned int i=0; i<polyctx.n; i++)
    printf("%lu ", mpi_to_u64(a.coeffs[i]));
  printf("\n");
  poly_rns_free(&ahat);
TEST_DONE();

  poly_mpi_free(&a);

TEST_DO("exit");
  polyctx_exit();
  mpi_release(Q);
TEST_DONE();

  mpi_release(q);
TEST_END();
}

int main()
{
  test_crt();
  return 0;
}
