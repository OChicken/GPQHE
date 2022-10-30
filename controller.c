/*
 * Control system parameters precomputation.
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

#include "controller.h"

#include <string.h>       /* memcpy, strcpy, strerror */

/* error */
#include <assert.h>
#include <errno.h>

/* serialization */
#include <fcntl.h>   /* open & creat, return int */
#include <unistd.h>  /* read & write, return ssize_t */

/* BLAS */
#include <cblas.h>
#include <lapacke.h>

void linear_controller_ctx_build(struct controller_ctx *controller)
{
  uint32_t n    = GPQHE_DIMX;
  uint32_t m    = GPQHE_DIMU;
  uint32_t N    = GPQHE_HORIZON;
  uint32_t Np1  = N+1;
  uint32_t nNp1 = n*(N+1);
  uint32_t mN   = m*N;
  int err;
  memcpy(controller->A, (double[]){
    0.77,-0.35,
    0.49, 0.91
  }, n*n*sizeof(double));
  memcpy(controller->B, (double[]){
    0.04,
    0.15
  }, n*m*sizeof(double));
  memcpy(controller->Q, (double[]){
    50, 25,
    25, 10
  }, n*n*sizeof(double));
  memcpy(controller->R, (double[]){
    1.0,
  }, m*m*sizeof(double));
  memcpy(controller->P, (double[]){
    150,  0,
      0, 10
  }, n*n*sizeof(double));

  double An[GPQHE_DIMX*GPQHE_DIMX] = {0}; /*  n*n    */
  double AA[(GPQHE_DIMX*(GPQHE_HORIZON+1)) * GPQHE_DIMX] = {0}; /* n(N+1)*n */
  double BB[(GPQHE_DIMX*(GPQHE_HORIZON+1)) * (GPQHE_DIMU*GPQHE_HORIZON)] = {0}; /* n(N+1)*mN */
  double QQ[(GPQHE_DIMX*(GPQHE_HORIZON+1)) * (GPQHE_DIMX*(GPQHE_HORIZON+1))] = {0}; /* n(N+1)*n(N+1) */
  double RR[(GPQHE_DIMU*GPQHE_HORIZON) * (GPQHE_DIMU*GPQHE_HORIZON)] = {0}; /* mN*mN */
  /* k=0 */
  for (uint32_t i=0; i<n*n; i+=n)
    An[i++]=1;
  for (uint32_t j=0, allA=n*n; j<allA; j++)
    AA[(j/n)*n+j%n]=An[j];
  /* k = 1 ~ N+1 */
  for (uint32_t k=1; k<Np1; k++)
  {
    double tmpBB[GPQHE_DIMX*GPQHE_DIMU] = {0}; /* n*m */
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
      n,m,n,1,An,n,controller->B,m,1,tmpBB,m);
    for (uint32_t i=k; i<Np1; i++)
    {
      for (uint32_t j=0, allB=n*m; j<allB; j++)
        BB[(i*n+j/m)*mN+(i-k)*m+j%m]=tmpBB[j];
    }
    double tmpAA[GPQHE_DIMX*GPQHE_DIMX] = {0}; /* n*n */
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
      n,n,n,1,An,n,controller->A,n,1,tmpAA,n);
    cblas_dcopy(n*n,tmpAA,1,An,1);
    for (uint32_t j=0, allA=n*n; j<allA; j++)
      AA[k*allA+j]=An[j];
  }
  for (uint32_t i=0; i<N; i++)
  {
    for (uint32_t j=0; j<n*n; j++)
      QQ[(i*n+j/n)*nNp1 + i*n+j%n]=controller->Q[j];
    for (uint32_t j=0; j<m*m; j++)
      RR[(i*m+j/m)*mN + i*m+j%m]=controller->R[j];
  }
  for (uint32_t j=0; j<n*n; j++)
    QQ[N*(n*nNp1)+N*n+(j/n)*nNp1+(j%n)]=controller->P[j];
  /* begin calculate K */
  double BTQ[(GPQHE_DIMU*GPQHE_HORIZON)*(GPQHE_DIMX*(GPQHE_HORIZON+1))] = {0}; /*  mN*n(N+1) */
  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,
    mN,nNp1,nNp1,1,BB,mN,QQ,nNp1,1,BTQ,nNp1);
  double BTQBR[(GPQHE_DIMU*GPQHE_HORIZON)*(GPQHE_DIMU*GPQHE_HORIZON)] = {0}; /* mN*mN */
  cblas_dcopy(mN*mN,RR,1,BTQBR,1);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
    mN,mN,nNp1,1,BTQ,nNp1,BB,mN,1,BTQBR,mN);
  double BTQBR_inv[(GPQHE_DIMU*GPQHE_HORIZON)*(GPQHE_DIMU*GPQHE_HORIZON)] = {0}; /* mN*mN */
  cblas_dcopy(mN*mN,BTQBR,1,BTQBR_inv,1);
  int ipiv[mN];
  err = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,mN,mN,BTQBR_inv,mN,ipiv);
  if(err)
    printf("%s\n",strerror(err));
  err = LAPACKE_dgetri(LAPACK_ROW_MAJOR,mN,BTQBR_inv,mN,ipiv);
  if(err)
    printf("%s\n",strerror(err));
  double BTQA[(GPQHE_DIMU*GPQHE_HORIZON)*GPQHE_DIMX] = {0}; /* mN*n */
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
    mN,n,nNp1,1,BTQ,nNp1,AA,n,1,BTQA,n);
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
    mN,n,mN,1,BTQBR_inv,mN,BTQA,n,1,controller->K,n);
  cblas_dscal(mN*n, -1, controller->K, 1);
#if 0
  /* print K */
  printf("K=\n");
  for (uint32_t i=0; i<m*N*n; i++)
  {
    if ((i%n==0)&(i!=0)) printf("\n");
    printf("% -12.8e ", controller->K[i]);
  }
  printf("\n");
#endif
}

void linear_controller_ctx_save(struct controller_ctx *controller)
{
  uint32_t n    = GPQHE_DIMX;
  uint32_t m    = GPQHE_DIMU;
  uint32_t N    = GPQHE_HORIZON;
  uint32_t mN   = m*N;

  int openfd;
  ssize_t writefd;
  openfd = open("precomp/A", O_WRONLY | O_CREAT, 0664);
  writefd = write(openfd, controller->A, n*n*sizeof(double));
  openfd = open("precomp/B", O_WRONLY | O_CREAT, 0664);
  writefd = write(openfd, controller->B, n*m*sizeof(double));
  openfd = open("precomp/K", O_WRONLY | O_CREAT, 0664);
  writefd = write(openfd, controller->K, mN*n*sizeof(double));
  close(writefd);
  close(openfd);
}

void linear_controller_ctx_load(struct controller_ctx *controller)
{
  uint32_t n    = GPQHE_DIMX;
  uint32_t m    = GPQHE_DIMU;
  uint32_t N    = GPQHE_HORIZON;
  uint32_t mN   = m*N;

  int openfd;
  ssize_t readfd;
  openfd = open(PATH_OF_PRECOMP "A", O_RDONLY, 0000);
  readfd = read(openfd, controller->A, n*n*sizeof(uint64_t));
  openfd = open(PATH_OF_PRECOMP "B", O_RDONLY, 0000);
  readfd = read(openfd, controller->B, n*m*sizeof(uint64_t));
  openfd = open(PATH_OF_PRECOMP "K", O_RDONLY, 0000);
  readfd = read(openfd, controller->K, mN*n*sizeof(uint64_t));
  close(readfd);
  close(openfd);
}
#if 0
void linear_controller_actuate(const char *name, double *x0, double *uopt,
  struct controller_ctx *controller)
{
  uint32_t n = GPQHE_DIMX;
  uint32_t m = GPQHE_DIMU;
  uint32_t N = GPQHE_HORIZON;
  int openfd;
  ssize_t readfd;
  double x[GPQHE_DIMX] = {0};
  double u[GPQHE_DIMU] = {0};
  /* load matrices */
  double A[GPQHE_DIMX*GPQHE_DIMX] = {0};
  openfd = open("precomp/A", O_RDONLY, 0000);
  readfd = read(openfd, A, n*n*sizeof(double));
  double B[GPQHE_DIMX*GPQHE_DIMU] = {0};
  openfd = open("precomp/B", O_RDONLY, 0000);
  readfd = read(openfd, B, n*m*sizeof(double));
  close(readfd);
  close(openfd);
  /* state evolution */
  double record[(GPQHE_DIMX+GPQHE_DIMU)*(GPQHE_HORIZON+1)] = {0};
  x[0]=x0[0],x[1]=x0[1];
  for (uint32_t i=0; i<N; i++)
  {
    for (uint32_t j=0; j<n; j++)
      record[i*(n+m)+j]=x[j];
    for (uint32_t j=0; j<m; j++)
      record[i*(n+m)+n+j]=u[j]=uopt[i*m+j];
    double Ax[GPQHE_DIMX] = {0};
    cblas_dgemv(CblasRowMajor,CblasNoTrans,n,n,1,A,n,x,1,1,Ax,1);
    double Bu[GPQHE_DIMX] = {0};
    cblas_dgemv(CblasRowMajor,CblasNoTrans,n,m,1,B,m,u,1,1,Bu,1);
    cblas_dcopy(n,Ax,1,x,1);
    cblas_daxpy(n,1,Bu,1,x,1);
  }
  for (uint32_t j=0; j<n; j++)
    record[N*(n+m)+j]=x[j];
  openfd = open(name, O_WRONLY | O_CREAT, 0664);
  ssize_t writefd = write(openfd, record, (n+m)*(N+1)*sizeof(double));
  close(writefd);
  close(openfd);
}
#endif