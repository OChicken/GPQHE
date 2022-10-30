/*
 * Number Theory Development Kit based on Libgcrypt
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

/* begin ntdk.c */

#include "ntdk.h"

/*******************************************************************************
 * Conversion                                                                  *
 ******************************************************************************/

void
show_mpi (MPI a)
{
  gcry_error_t err = GPG_ERR_NO_ERROR;
  gcry_sexp_t data;
  char *buf;
  size_t size;
  err = gcry_sexp_build(&data, NULL, "%M", a);
  if (err)
    fprintf(stderr, "Error in %s.", __func__);
  size = gcry_sexp_sprint (data, GCRYSEXP_FMT_ADVANCED, NULL, 0);
  buf = (char *)malloc (size);
  gcry_sexp_sprint (data, GCRYSEXP_FMT_ADVANCED, buf, size);
  fflush(stdout);
  fprintf (stderr, "%s", buf);
  free (buf);
  gcry_sexp_release(data);
}

char *
sprint_hex(MPI N)
{
  gcry_error_t err = GPG_ERR_NO_ERROR;
  unsigned char *pc;
  char buf[BUFSIZ];
  size_t buflen;
  err = gcry_mpi_aprint(GCRYMPI_FMT_USG, &pc, &buflen, N);
  if (err)
    fprintf(stderr, "Error in %s.", __func__);
  char *p=(char *)malloc(buflen*sizeof(char *)+1);
  strcpy(p,"");
  if (gcry_mpi_cmp_ui(N, 0)<0)
    sprintf(p, "-");
  else if (gcry_mpi_cmp_ui(N, 0)==0)
    sprintf(p, "0");
  for (unsigned char *s=pc; buflen; buflen--, s++)
  {
    sprintf(buf, "%02x", *s);
    strcat(p, buf);
  }
  return p;
}

int64_t
gcry_mpi_to_int64 (MPI a)
{
  int64_t num;
  size_t length = gcry_mpi_get_nbits (a);
  /* if is zero */
  if (!length)
    return 0;
  /* if not zero, and neg */
  if (gcry_mpi_is_neg(a))
  {
    MPI b = gcry_mpi_new(0);
    gcry_mpi_set(b,a);
    gcry_mpi_neg(b,b);
    num = gcry_mpi_test_bit(b, length);
    while (length-- > 0)
    {
      num <<= 1;
      num += gcry_mpi_test_bit(b, length);
    }
    num = -num;
    gcry_mpi_release(b);
  }
  /* if not zero, and pos */
  else
  {
    num = gcry_mpi_test_bit(a, length);
    while (length-- > 0)
    {
      num <<= 1;
      num += gcry_mpi_test_bit(a, length);
    }
  }
  return num;
}

uint64_t
gcry_mpi_to_uint64 (MPI a)
{
  uint64_t num;
  size_t length = gcry_mpi_get_nbits (a);
  /* if is zero */
  if (!length)
    return 0;
  /* if not zero, and pos */
  assert(gcry_mpi_is_neg(a)==0);
  num = gcry_mpi_test_bit(a, length);
  while (length-- > 0)
  {
    num <<= 1;
    num += gcry_mpi_test_bit(a, length);
  }
  return num;
}

__int128_t
gcry_mpi_to_int128 (MPI a)
{
  __int128_t num;
  size_t length = gcry_mpi_get_nbits (a);
  /* if is zero */
  if (!length)
    return 0;
  /* if not zero, and neg */
  if (gcry_mpi_is_neg(a))
  {
    MPI b = gcry_mpi_new(0);
    gcry_mpi_set(b,a);
    gcry_mpi_neg(b,b);
    num = gcry_mpi_test_bit(b, length);
    while (length-- > 0)
    {
      num <<= 1;
      num += gcry_mpi_test_bit(b, length);
    }
    num = -num;
    gcry_mpi_release(b);
  }
  /* if not zero, and pos */
  else
  {
    num = gcry_mpi_test_bit(a, length);
    while (length-- > 0)
    {
      num <<= 1;
      num += gcry_mpi_test_bit(a, length);
    }
  }
  return num;
}

__uint128_t
gcry_mpi_to_uint128 (MPI a)
{
  __uint128_t num;
  size_t length = gcry_mpi_get_nbits (a);
  /* if is zero */
  if (!length)
    return 0;
  /* if not zero, and pos */
  assert(gcry_mpi_is_neg(a)==0);
  num = gcry_mpi_test_bit(a, length);
  while (length-- > 0)
  {
    num <<= 1;
    num += gcry_mpi_test_bit(a, length);
  }
  return num;
}

MPI
gcry_int64_to_mpi (int64_t a)
{
  MPI num = gcry_mpi_new(0);
  if (a<0)
  {
    gcry_mpi_set_ui(num, -a);
    gcry_mpi_neg(num, num);
  }
  else
    gcry_mpi_set_ui(num, a);
  return num;
}

MPI
gcry_int128_to_mpi (__int128_t a)
{
  MPI num = gcry_mpi_new(0);
  if (a<0)
  {
    gcry_mpi_set_ui(num, -a);
    gcry_mpi_neg(num, num);
  }
  else
    gcry_mpi_set_ui(num, a);
  return num;
}

uint64_t mpi_buf2uint64 (const char *buf, const size_t len)
{
  assert(len<=8*sizeof(uint64_t));
  uint64_t num = 0;
  for (size_t i=0; i<len; i++)
    num += ((int8_t)buf[i]<0? 0:1)<<(len-1-i);
  return num;
}

uint64_t
ntdk_buf_to_u64 (const unsigned char *buf, const size_t l)
{
  assert(l<=64);
  uint64_t num=0;
  for (size_t k=0; k<l; k++)
  {
    size_t i=k/8;
    size_t j=k%8;
    num += ((uint64_t)((buf[i]>>j)&0x01))<<k;
  }
  return num;
}

/*******************************************************************************
 * Arithmetic Utils                                                            *
 ******************************************************************************/

MPI
gcry_mpi_sqrt(MPI x)
{
  /* 程序员面试题精解（2）— 平方根运算 | 网络热度
   * https://www.packetmania.net/2021/07/23/PGITVW-2-sqrt/ */
  assert(gcry_mpi_is_neg(x)==0);
  MPI low = gcry_mpi_new(0);
  MPI high = gcry_mpi_new(0);
  MPI mid = gcry_mpi_new(0);
  MPI tmp = gcry_mpi_new(0);
  MPI two = gcry_mpi_set_ui(NULL, 2);
  if (gcry_mpi_cmp_ui(x, 1)>0)
  {
    gcry_mpi_set_ui(low, 1);
    gcry_mpi_set(high, x);
  }
  else
  {
    gcry_mpi_set(low, x);
    gcry_mpi_set_ui(high, 1);
  }
  while (gcry_mpi_cmp(high, low)>=0)
  {
    gcry_mpi_add(mid, high, low);
    gcry_mpi_div(mid, NULL, mid, two, -1); /* mid=(high+low)/2 */
    if (gcry_mpi_cmp(mid, low)==0)
      break;
    gcry_mpi_mul(tmp, mid, mid);
    if (gcry_mpi_cmp(tmp,x)>0)
      gcry_mpi_set(high, mid);// high = mid;
    else
      gcry_mpi_set(low, mid);//low = mid;
  }
  gcry_mpi_release(tmp);
  gcry_mpi_release(low);
  gcry_mpi_release(high);
  return mid;
}

size_t
ntdk_log2 (unsigned long x)
{
  size_t w = 0;
  while (x>1)
  {
    w += 1;
    x >>= 1;
  }
  return w;
}

size_t
ntdk_get_nbits(uint64_t x)
{
  return ntdk_log2(x)+1;
}

long
mpi_br_width (long x, size_t w)
{
  long y = 0;
  while (w>0){
    if ( (x&1) == 1)
      y += 1;
    if (w>1)
      y <<= 1;
    x >>= 1;
    w -= 1;
  }
  return y;
}

uint32_t
ntdk_br_u32(uint32_t x)
{
  x = ((x & 0xaaaaaaaa)>> 1) | ((x & 0x55555555)<< 1);
  x = ((x & 0xcccccccc)>> 2) | ((x & 0x33333333)<< 2);
  x = ((x & 0xf0f0f0f0)>> 4) | ((x & 0x0f0f0f0f)<< 4);
  x = ((x & 0xff00ff00)>> 8) | ((x & 0x00ff00ff)<< 8);
  x = ((x & 0xffff0000)>>16) | ((x & 0x0000ffff)<<16);
  return x;
}

uint64_t
ntdk_br_u64(uint64_t x)
{
  x = ((x & 0xaaaaaaaaaaaaaaaa)>> 1) | ((x & 0x5555555555555555)<< 1);
  x = ((x & 0xcccccccccccccccc)>> 2) | ((x & 0x3333333333333333)<< 2);
  x = ((x & 0xf0f0f0f0f0f0f0f0)>> 4) | ((x & 0x0f0f0f0f0f0f0f0f)<< 4);
  x = ((x & 0xff00ff00ff00ff00)>> 8) | ((x & 0x00ff00ff00ff00ff)<< 8);
  x = ((x & 0xffff0000ffff0000)>>16) | ((x & 0x0000ffff0000ffff)<<16);
  x = ((x & 0xffffffff00000000)>>32) | ((x & 0x00000000ffffffff)<<32);
  return x;
}

uint64_t
ntdk_sqrt_u64(uint64_t x)
{
  uint64_t low, high, mid, tmp;
  if (x>1)
  {
    low = 1;
    high = x;
  }
  else
  {
    low = x;
    high = 1;
  }
  while (high>=low)
  {
    mid = (high+low)/2;
    if (mid==low)
      break;
    tmp = mid*mid;
    if (tmp>x)
      high = mid;
    else
      low = mid;
  }
  return mid;
}

uint64_t
ntdk_pow(uint64_t b, uint64_t e)
{
  uint64_t ret = 1;
  while (e>0)
  {
    if (e&1)
      ret *= b;
    e >>= 1;
    b *= b;
  }
  return ret;
}

uint64_t
ntdk_inv_u64(uint64_t x)
{
  return ntdk_pow(x, (uint64_t)(-1));
}

uint64_t
ntdk_mulm(uint64_t a, uint64_t b, uint64_t m)
{
  __uint128_t mul = (__uint128_t)(a) * (__uint128_t)(b);
  mul %= (__uint128_t)(m);
  return (uint64_t)(mul);
}

uint64_t
ntdk_mulm_barrett(uint64_t a, uint64_t b, uint64_t p, uint64_t pr)
{
  __uint128_t mul = (__uint128_t)(a) * (__uint128_t)(b);
  uint64_t abot = (uint64_t)(mul);
  uint64_t atop = (uint64_t)(mul>>64);
  __uint128_t tmp;
  tmp   = (__uint128_t)(abot) * (__uint128_t)(pr);
  tmp >>= 64;
  tmp  += (__uint128_t)(atop) * (__uint128_t)(pr);
  tmp >>= 2*(PBND+1) - 64;
  tmp  *= p;
  tmp   = mul - tmp;
  uint64_t r = (uint64_t)(tmp);
  if (r>=p) r -= p;
  return r;
}

uint64_t
ntdk_powm(uint64_t b, uint64_t e, uint64_t m)
{
  uint64_t ret = 1;
  while (e>0)
  {
    if (e&1)
      ret = ntdk_mulm(ret,b,m);
    e >>= 1;
    b = ntdk_mulm(b,b,m);
  }
  return ret;
}

void
ntdk_smod(MPI *a, const MPI m)
{
  MPI zero = gcry_mpi_set_ui(NULL, 0);
  MPI m_half = gcry_mpi_new(0);
  gcry_mpi_div(m_half, NULL, m, gcry_mpi_set_ui(NULL, 2), -1);
  while (gcry_mpi_cmp(*a, zero)<0)
    gcry_mpi_add(*a, *a, m);
  gcry_mpi_mod(*a, *a, m);
  if (gcry_mpi_cmp(*a, m_half)>0)
    gcry_mpi_sub(*a, *a, m);
}

void
ntdk_div_round(MPI q, MPI dividend, MPI divisor)
{
  assert(!gcry_mpi_is_neg(divisor));
  uint8_t neg=0;
  if (gcry_mpi_is_neg(dividend))
    neg=1;
  MPI divisor_half = gcry_mpi_new(0);
  MPI r = gcry_mpi_new(0);
  gcry_mpi_div(divisor_half, NULL, divisor, GCRYMPI_CONST_TWO, 0);
  gcry_mpi_div(q, r, dividend, divisor, -1);
  if ((!neg) && (gcry_mpi_cmp(r, divisor_half)>0))
    gcry_mpi_add(q, q, GCRYMPI_CONST_ONE);
  if (neg && (gcry_mpi_cmp_ui(r,0)!=0))
  {
    //gcry_mpi_neg(q, q);
    if (gcry_mpi_cmp(r, divisor_half)>0)
      gcry_mpi_add(q,q,GCRYMPI_CONST_ONE);
  }
}

uint64_t
ntdk_m_root_of_unity(uint64_t m, uint64_t p)
{
  uint64_t phi = p-1;
  /* factorize */
  uint64_t N = phi;
  uint64_t *factors = (uint64_t *)calloc(
    ntdk_get_nbits(ntdk_sqrt_u64(phi))+1, sizeof(uint64_t));
  size_t i=0;
  while (N%2==0)
  {
    factors[i++] = 2;
    N /= 2;
  }
  for (uint64_t n=3; n<ntdk_sqrt_u64(N); n++)
  {
    while (N%n==0)
    {
      factors[i++] = n;
      N /= n;
    }
  }
  if (N>2)
    factors[i] = N;
  /* generator */
  uint64_t g;
  for (g=2; g<=phi; g++)
  {
    uint8_t flag = 0;
    for (size_t j=0; factors[j]; j++)
    {
      if (ntdk_powm(g, phi/factors[j], p)==1)
      {
        flag = 1;
        break;
      }
    }
    if (flag==0)
      break;
  }
  free(factors);
  /* m-th root of unity */
  assert(phi%m==0);
  return ntdk_powm(g, phi/m, p);
}



/*******************************************************************************
 * Polynomial                                                                  *
 * ****************************************************************************/

void
poly_ctx_build(poly_ctx *poly, size_t d)
{
  poly->d = d;
  poly->m = 2*d;
  uint64_t logd = ntdk_log2(d);
  size_t np = (2 + logd + 4*LOGQ + PBND - 1) / PBND;
  //MPI two    = gcry_mpi_set_ui(NULL,2);
  poly->p_prod      = gcry_calloc(np, sizeof(MPI));
  poly->ph_prod     = gcry_calloc(np, sizeof(MPI));
  poly->ph          = gcry_calloc(np, sizeof(MPI));
  poly->ph_inv_modp = gcry_calloc(np*np, sizeof(MPI));

  /* set primes */
  uint64_t prime = (1ull<<PBND) + 1;
  for (size_t i = 0; i < np; i++)
  {
    while(1)
    {
      prime += poly->m;
      if(!gcry_prime_check(gcry_mpi_set_ui(NULL,prime), 0))
      {
        poly->p[i] = prime;
        break;
      }
    }
    poly->p_prod[i]      = gcry_mpi_new(0);
    poly->ph_prod[i]     = gcry_mpi_new(0);
    poly->ph[i]          = gcry_mpi_new(0);
    poly->ph_inv_modp[i] = gcry_mpi_new(0);
    if (i==0)
      gcry_mpi_set_ui(NULL, poly->p[i]);
    else
      gcry_mpi_mul(poly->p_prod[i], poly->p_prod[i-1], gcry_mpi_set_ui(NULL, poly->p[i]));
  }

  poly->d_inv_scaled       = (uint64_t *)calloc(d, sizeof(uint64_t));
  poly->rootpow_scaled     = (uint64_t **)calloc(d, sizeof(uint64_t *));
  poly->rootpow_inv_scaled = (uint64_t **)calloc(d, sizeof(uint64_t *));
  for (size_t i=0; i<np; i++)
  {
    //red_ss_array[i] = _ntl_general_rem_one_struct_build(pVec[i]);
    uint64_t pi = poly->p[i];
    poly->p_inv[i] = ntdk_inv_u64(pi);
    poly->pr[i]    = ((__uint128_t)(1) << (2*(PBND+1))) / pi;
    uint64_t root     = ntdk_m_root_of_unity(poly->m, pi);
    uint64_t root_inv = ntdk_powm(root, pi-2, pi);
    uint64_t d_inv    = ntdk_powm(poly->d, pi-2, pi);
    poly->d_inv_scaled[i]       = ntdk_mulm(ntdk_mulm(d_inv, 1ull<<32, pi), 1ull<<32, pi);
    poly->rootpow_scaled[i]     = (uint64_t *)calloc(d, sizeof(uint64_t));
    poly->rootpow_inv_scaled[i] = (uint64_t *)calloc(d, sizeof(uint64_t));
    uint64_t power     = 1;
    uint64_t power_inv = 1;
    for (size_t j=0; j<d; j++) {
      uint32_t pj = ntdk_br_u32((uint32_t)(j)) >> (32-logd);
      uint64_t rootpow     = power;
      uint64_t rootpow_inv = power_inv;
      poly->rootpow_scaled[i][pj]     = ntdk_mulm(ntdk_mulm(rootpow, 1ull<<32, pi), 1ull<<32, pi);
      poly->rootpow_inv_scaled[i][pj] = ntdk_mulm(ntdk_mulm(rootpow_inv, 1ull<<32, pi), 1ull<<32, pi);
      power     = ntdk_mulm(power    , root, pi);
      power_inv = ntdk_mulm(power_inv, root_inv, pi);
    }
  }
}

/*******************************************************************************
 * BLAS                                                                        *
 ******************************************************************************/

void
blas_lgemv(size_t m, size_t n, long *y, long *A, long *x)
{
  for (size_t i=0; i<m; i++)
  {
    y[i] = 0;
    for (size_t j=0; j<n; j++)
      y[i] += A[i*n+j]*x[j];
  }
}

void
blas_zgemv(size_t m, size_t n, _Complex double *y, _Complex double *A, _Complex double *x)
{
  for (size_t i=0; i<m; i++)
  {
    y[i] = 0;
    for (size_t j=0; j<n; j++)
      y[i] += A[i*n+j]*x[j];
  }
}

void
blas_xpy(size_t m, long *r, long *x, long *y)
{
  for (size_t i=0; i<m; i++)
    r[i] = x[i]+y[i];
}

void
blas_zxpy(size_t m, _Complex double *r, _Complex double *x, _Complex double *y)
{
  for (size_t i=0; i<m; i++)
    r[i] = x[i]+y[i];
}

void
blas_zxsy(size_t m, _Complex double *r, _Complex double *x, _Complex double *y)
{
  for (size_t i=0; i<m; i++)
    r[i] = x[i]-y[i];
}

void
blas_xmy(size_t m, long *r, long *x, long *y)
{
  for (size_t i=0; i<m; i++)
    r[i] = x[i]*y[i];
}

void
blas_zxmy(size_t m, _Complex double *r, _Complex double *x, _Complex double *y)
{
  for (size_t i=0; i<m; i++)
    r[i] = x[i]*y[i];
}

void
blas_scal(size_t m, long *y, long alpha, long *x)
{
  for (size_t i=0;i<m; i++)
    y[i] = alpha*x[i];
}

void
blas_ldiag(size_t m, long *diag, size_t index, long *A)
{
  for (size_t i=0; i<m; i++)
    diag[i] = A[(i%m)*m + (index+i)%m];
}

void
blas_ddiag(size_t m, double *diag, size_t index, double *A)
{
  for (size_t i=0; i<m; i++)
    diag[i] = A[(i%m)*m + (index+i)%m];
}

void
blas_zdiag(size_t m, _Complex double *diag, size_t index, _Complex double *A)
{
  for (size_t i=0; i<m; i++)
    diag[i] = A[(i%m)*m + (index+i)%m];
}

void
blas_lrot(size_t m, long *y, int64_t rot, long *x)
{
  assert(y!=x);
  for (size_t i=0; i<m; i++)
  {
    int64_t idx = (i+rot)%m;
    if (idx<0)
      idx += m;
    assert(idx>=0);
    y[i] = x[idx];
  }
}

void
blas_drot(size_t m, double *y, int64_t rot, double *x)
{
  assert(y!=x);
  for (size_t i=0; i<m; i++)
  {
    int64_t idx = (i+rot)%m;
    if (idx<0)
      idx += m;
    assert(idx>=0);
    y[i] = x[idx];
  }
}

void
blas_zrot(size_t m, _Complex double *y, int64_t rot, _Complex double *x)
{
  assert(y!=x);
  for (size_t i=0; i<m; i++)
  {
    int64_t idx = (i+rot)%m;
    if (idx<0)
      idx += m;
    assert(idx>=0);
    y[i] = x[idx];
  }
}

void
blas_conj(size_t m, _Complex double *y, _Complex double *x)
{
  for (size_t i=0; i<m; i++)
    y[i] = conj(x[i]);
}

void
blas_trans(size_t m, size_t n, long *AT, long *A)
{
  for (size_t i=0; i<m; i++)
  {
    for (size_t j=0; j<n; j++)
      AT[j*m+i] = A[i*n+j];
  }
}

uint8_t
blas_mpivec_equal(size_t m, MPI *y, MPI *x)
{
  for (size_t i=0; i<m; i++)
  {
    if (gcry_mpi_cmp(y[i],x[i])!=0)
      return 0;
  }
  return 1;
}

uint8_t
blas_lvec_equal(size_t m, long *y, long *x)
{
  for (size_t i=0; i<m; i++)
  {
    if (y[i]!=x[i])
      return 0;
  }
  return 1;
}

uint8_t
blas_zvec_equal(size_t m, _Complex double *y, _Complex double *x, double threshold)
{
  for (size_t i=0; i<m; i++)
  {
    double diff = y[i]-x[i];
    diff = (diff>=0? diff : -diff);
    if (diff>threshold)
      return 0;
  }
  return 1;
}

double
blas_dznrmmax(size_t m, _Complex double *z)
{
  double tmp = 0;
  double norm = 0;
  for (size_t i=0; i<m; i++)
  {
    tmp = cabs(z[i]);
    if (norm<tmp)
      norm = tmp;
  }
  return norm;
}

double
blas_dznrmmax_dist(size_t m, _Complex double *y, _Complex double *x)
{
  double tmp = 0;
  double dist = 0;
  for (size_t i=0; i<m; i++) {
    tmp = cabs(y[i] - x[i]);
    if (dist<tmp)
      dist = tmp;
  }
  return dist;
}

double blas_dnrmmax_dist(size_t m, double *y, double *x)
{
  double tmp = 0;
  double dist = 0;
  for (size_t i=0; i<m; i++) {
    tmp = fabs(y[i] - x[i]);
    if (dist<tmp)
      dist = tmp;
  }
  return dist;
}

/*******************************************************************************
 * Sample                                                                      *
 ******************************************************************************/

void
sample_z01vec(size_t m, _Complex double *vec)
{
  unsigned char buf[m*2];
  gcry_randomize(buf, 2*m, GCRY_STRONG_RANDOM);
  for (size_t i=0; i<m; i++)
    vec[i] = (double)buf[i]/256 + I*(double)buf[i+m]/256;
}

/*
void
sample_cbd2_eta1(size_t m, int16_t *vec)
{
  unsigned char buf[KYBER_ETA1*KYBER_N/4];
  gcry_randomize(buf, sizeof(buf), GCRY_WEAK_RANDOM);
  unsigned int i,j;
  uint32_t t,d;
  int16_t a,b;

  for(i=0;i<KYBER_N/8;i++) {
    t  = load32_littleendian(buf+4*i);
    d  = t & 0x55555555;
    d += (t>>1) & 0x55555555;

    for(j=0;j<8;j++) {
      a = (d >> (4*j+0)) & 0x3;
      b = (d >> (4*j+2)) & 0x3;
      r->coeffs[8*i+j] = a - b;
    }
  }
}
*/

void
sample_discrete_gaussian(size_t m, int16_t *vec, double sigma)
{
  for (size_t i=0; i<m; i+=2)
  {
    unsigned char buf[2];
    gcry_randomize(buf, 2, GCRY_STRONG_RANDOM);
    double r1 = (double)buf[0]/256;
    double r2 = (double)buf[1]/256;
    double theta = 2*M_PI*r1;
    double rr = sqrt(-2*log(r2)) * sigma;
    vec[i  ] = (int16_t)floor(rr*cos(theta)+0.5);
    vec[i+1] = (int16_t)floor(rr*sin(theta)+0.5);
  }
}

void
sample_hwt(size_t m, int8_t *vec, size_t h)
{
  MPI tmp = gcry_mpi_new(0);
  do {
    gcry_mpi_randomize(tmp, h, GCRY_STRONG_RANDOM);
  } while(gcry_mpi_get_nbits(tmp)!=h);
  unsigned int logm = ntdk_log2(m);
  size_t idx = 0; /* total nonzero entries count */
  while (idx<h)
  {
    char buf[logm];
    gcry_randomize(buf, logm, GCRY_WEAK_RANDOM);
    int64_t i = mpi_buf2uint64(buf, logm);
    if (vec[i]==0)
    {
      vec[i] = (gcry_mpi_test_bit(tmp, idx)==0)? 1 : -1;
      idx++;
    }
  }
  gcry_mpi_release(tmp);
}

void
sample_zero_center(size_t m, int16_t *vec)
{
  MPI tmp = gcry_mpi_new(0);
  do
  {
    gcry_mpi_randomize(tmp, 2*m, GCRY_STRONG_RANDOM);
  } while(gcry_mpi_get_nbits(tmp)!=2*m);
  for (size_t i=0; i<m; i++)
    vec[i] = (gcry_mpi_test_bit(tmp, 2*i)==0)? 0 : (gcry_mpi_test_bit(tmp, 2*i+1)==0)? 1: -1;
  gcry_mpi_release(tmp);
}

void
sample_uniform_u64(size_t m, uint64_t *r, const uint64_t q)
{
  size_t l = ntdk_get_nbits(q);
  size_t buflen = l/8+1;
  unsigned char *buf = (unsigned char *)malloc(buflen*sizeof(unsigned char));
  for (size_t i=0; i<m; i++)
  {
    gcry_randomize(buf, buflen, GCRY_STRONG_RANDOM);
    r[i] = ntdk_buf_to_u64(buf,l);
    while (r[i]>=q)
      r[i]-=q;
  }
  free(buf);
}




/*******************************************************************************
 * CRT                                                                         *
 ******************************************************************************/

void
crt_build(crt_ctx *crt, size_t size_p, size_t num_p)
{
  crt->L = num_p;
  gcry_error_t err = GPG_ERR_NO_ERROR;
  crt->p = gcry_calloc(crt->L, sizeof(MPI));
  crt->g = gcry_calloc(crt->L, sizeof(MPI));
  crt->m = gcry_calloc(crt->L, sizeof(MPI));
  crt->c = gcry_calloc(crt->L, sizeof(MPI));
  crt->M = gcry_mpi_set_ui(NULL, 1);
  for (size_t i=0; i<crt->L; i++)
  {
    crt->p[i] = gcry_mpi_new(0);
    crt->g[i] = gcry_mpi_new(0);
    crt->m[i] = gcry_mpi_new(0);
    crt->c[i] = gcry_mpi_new(0);
  }
  MPI mpi_p = gcry_mpi_new(0);
  /* generate primes and generators */
  MPI *factors = NULL;
  for (size_t i=0; i<crt->L; i++)
  {
    err = gcry_prime_generate(&crt->p[i], size_p,
          0 /* factor_bits */, &factors /* **factors */,
          NULL, NULL, /* cb_func, cb_arg */
          GCRY_WEAK_RANDOM, /* random_level */
          0 /* flags */);
    if (err || (gcry_mpi_get_nbits(crt->p[i])<size_p))
      fprintf(stderr, "Error in %s: generate primes fail.", __func__);
    err = gcry_prime_group_generator (&crt->g[i], crt->p[i], factors, NULL);
    if (err)
      fprintf(stderr, "Error in %s: generate primes fail.", __func__);
    gcry_mpi_mul(crt->M, crt->M, crt->p[i]);
    #if DEBUG == 1
    printf("g(%lu)=%lu\n", gcry_mpi_to_uint64(crt->p[i]), gcry_mpi_to_uint64(crt->g[i]));
    #endif
  }
  gcry_prime_release_factors (factors);
  factors = NULL;
  /* get mi and ci */
  MPI mi_inv = gcry_mpi_new(0);
  for (size_t i=0; i<crt->L; i++)
  {
    gcry_mpi_div(crt->m[i], NULL, crt->M, crt->p[i], -1); /* mi = M/p*/
    gcry_mpi_invm(mi_inv, crt->m[i], crt->p[i]);          /* inv(mi, pi) */
    gcry_mpi_mul(crt->c[i], crt->m[i], mi_inv);       /* ci = mi*inv(mi,pi) */
  }
  gcry_mpi_release(mi_inv);
  gcry_mpi_release(mpi_p);
}

void
crt_fwd(MPI *a, MPI n, crt_ctx *crt)
{
  for (size_t i=0; i<crt->L; i++)
  {
    gcry_mpi_mul(a[i], n, gcry_mpi_set_ui(NULL, 1));
    gcry_mpi_mod(a[i], a[i], crt->p[i]);
  }
}

void
crt_inv(MPI *n, MPI *a, crt_ctx *crt)
{
  MPI aici     = gcry_mpi_new(0);
  MPI aici_sum = gcry_mpi_set_ui(NULL,0);
  for (size_t i=0; i<crt->L; i++)
  {
    /* aici = ai*ci = ai*(mi*inv(mi,pi)) */
    gcry_mpi_mul(aici, a[i], crt->c[i]);
    gcry_mpi_add(aici_sum, aici_sum, aici);
  }
  gcry_mpi_mod(aici_sum, aici_sum, crt->M); /* aici_sum % M*/
  *n = gcry_mpi_copy(aici_sum);
  gcry_mpi_release(aici);
  gcry_mpi_release(aici_sum);
}

/*******************************************************************************
 * NTT                                                                         *
 ******************************************************************************/

void
ntt_build(ntt_ctx *ntt, size_t d, MPI p, MPI g)
{
  assert((d&(d-1))==0);
  assert(gcry_prime_check(p, 0)==0);
  ntt->d = d;
  ntt->p = gcry_mpi_copy(p);
  ntt->g = gcry_mpi_copy(g);
  ntt->zetas     = gcry_calloc(ntt->d, sizeof(MPI));
  ntt->zetas_inv = gcry_calloc(ntt->d, sizeof(MPI));
  ntt->zetas[0]     = gcry_mpi_set_ui(NULL, 1);
  ntt->zetas_inv[0] = gcry_mpi_set_ui(NULL, 1);
  /* M-th primitive root zeta & zeta_inv (M is cyclotomic index) */
  MPI phi = gcry_mpi_copy(p);
  gcry_mpi_sub_ui(phi, p, 1); /* phi = p-1 since modulus p is prime */
  MPI tmp      = gcry_mpi_new(0);
  MPI zeta     = gcry_mpi_new(0);
  MPI zeta_inv = gcry_mpi_new(0);
  gcry_mpi_div(tmp, NULL, phi, gcry_mpi_set_ui(NULL, 2*d), -1);
  gcry_mpi_powm(zeta, ntt->g, tmp, ntt->p); /* zeta = g^(phi/(2d)) mod p */
  gcry_mpi_invm(zeta_inv, zeta, ntt->p);
  /* powers of each primitive roots until degree d */
  for (size_t i=1; i<ntt->d; i++)
  {
    ntt->zetas[i] = gcry_mpi_new(0);
    gcry_mpi_mul(tmp, ntt->zetas[i-1], zeta);
    gcry_mpi_mod(ntt->zetas[i], tmp, p);
    ntt->zetas_inv[i] = gcry_mpi_new(0);
    gcry_mpi_mul(tmp, ntt->zetas_inv[i-1], zeta_inv);
    gcry_mpi_mod(ntt->zetas_inv[i], tmp, p);
  }
  gcry_mpi_release(phi);
  gcry_mpi_release(tmp);
  gcry_mpi_release(zeta);
  gcry_mpi_release(zeta_inv);
}

void
ntt_fwd(size_t m, MPI *X, const MPI *x, const MPI p, const ntt_ctx *ntt)
{
  assert(m==ntt->d);
  MPI *tmp = gcry_calloc(m, sizeof(MPI));
  MPI omega = gcry_mpi_new(0);
  MPI butterfly_plus  = gcry_mpi_new(0);
  MPI butterfly_minus = gcry_mpi_new(0);
  for (size_t i=0; i<m; i++)
  {
    tmp[i] = gcry_mpi_new(0);
    gcry_mpi_mul(tmp[i], x[i], ntt->zetas[i]);
    gcry_mpi_mod(tmp[i], tmp[i], p);
  }
  for (size_t i=0; i<m; i++)
    X[i] = tmp[mpi_br_width(i, ntdk_log2(m))];
  size_t width = ntdk_get_nbits(m);
  for (size_t logm=1; logm<width; logm++)
  {
    for (size_t j=0; j<m; j+=(1<<logm))
    {
      size_t half = 1<<(logm-1);
      for (size_t i=0; i<half; i++)
      {
        size_t idx_even = j+i;
        size_t idx_odd  = j+i+half;
        size_t zeta_idx = i<<(width-logm);
        gcry_mpi_mul(omega, ntt->zetas[zeta_idx], X[idx_odd]);
        gcry_mpi_mod(omega, omega, p);
        gcry_mpi_add(butterfly_plus, X[idx_even], omega);
        gcry_mpi_mod(butterfly_plus, butterfly_plus, p);
        gcry_mpi_sub(butterfly_minus, X[idx_even], omega);
        gcry_mpi_mod(butterfly_minus, butterfly_minus, p);
        X[idx_even] = gcry_mpi_copy(butterfly_plus);
        X[idx_odd]  = gcry_mpi_copy(butterfly_minus);
      }
    }
  }
  gcry_free(tmp);
  gcry_mpi_release(omega);
  gcry_mpi_release(butterfly_plus);
  gcry_mpi_release(butterfly_minus);
}

void
ntt_inv(size_t m, MPI *x, const MPI *X, const MPI p, const ntt_ctx *ntt)
{
  assert(m==ntt->d);
  MPI *tmp = gcry_calloc(m, sizeof(MPI));
  MPI omega = gcry_mpi_new(0);
  MPI butterfly_plus  = gcry_mpi_new(0);
  MPI butterfly_minus = gcry_mpi_new(0);
  for (size_t i=0; i<m; i++)
  {
    tmp[i] = gcry_mpi_new(0);
    tmp[i] = X[mpi_br_width(i, ntdk_log2(m))];
  }
  size_t width = ntdk_get_nbits(m);
  for (size_t logm=1; logm<width; logm++)
  {
    for (size_t j=0; j<m; j+=(1<<logm))
    {
      size_t half = 1<<(logm-1);
      for (size_t i=0; i<half; i++)
      {
        size_t idx_even = j+i;
        size_t idx_odd  = j+i+half;
        size_t zeta_idx = i<<(width-logm);
        gcry_mpi_mul(omega, ntt->zetas_inv[zeta_idx], tmp[idx_odd]);
        gcry_mpi_mod(omega, omega, p);
        gcry_mpi_add(butterfly_plus, tmp[idx_even], omega);
        gcry_mpi_mod(butterfly_plus, butterfly_plus, p);
        gcry_mpi_sub(butterfly_minus, tmp[idx_even], omega);
        gcry_mpi_mod(butterfly_minus, butterfly_minus, p);
        tmp[idx_even] = gcry_mpi_copy(butterfly_plus);
        tmp[idx_odd]  = gcry_mpi_copy(butterfly_minus);
      }
    }
  }
  MPI mpi_poly_d_inv = gcry_mpi_new(0);
  gcry_mpi_invm(mpi_poly_d_inv, gcry_mpi_set_ui(NULL, ntt->d), p);
  for (size_t i=0; i<m; i++)
  {
    gcry_mpi_mul(tmp[i], tmp[i], ntt->zetas_inv[i]);
    gcry_mpi_mul(tmp[i], tmp[i], mpi_poly_d_inv);
    gcry_mpi_mod(x[i], tmp[i], p);
  }
  gcry_free(tmp);
  gcry_mpi_release(omega);
  gcry_mpi_release(butterfly_plus);
  gcry_mpi_release(butterfly_minus);
  gcry_mpi_release(mpi_poly_d_inv);
}

void
fft_build(fft_ctx *fft, size_t d)
{
  /* The degree d polynomial is of complex coefficients, so there are altogether
     2*d real coefficients. Note that 2*d here is NOT cyclotomic index! */
  fft->L = 2*d;
  /* primitive root or order L */
  fft->zetas = (_Complex double *)calloc(fft->L, sizeof(_Complex double));
  /* powers of each primitive roots until the number of slots. */
  for (size_t i=0; i<fft->L; i++)
  {
    double theta = 2*M_PI*i / fft->L;
    fft->zetas[i] = cos(theta) + I*sin(theta);
  }
  /* rotation group */
  size_t num_slots = d/2;
  fft->rot_group = (int64_t *)calloc(num_slots, sizeof(int64_t));
  fft->rot_group[0] = 1;
  for (size_t i=1; i<num_slots; i++)
    fft->rot_group[i] = (5*fft->rot_group[i-1])%(fft->L);
}

void
fft_fwd(size_t m, _Complex double *X, const _Complex double *x, const fft_ctx *fft)
{
  /* m is not the length of fft_ctx->L. It should satisfy m <= 2*d = L.
     Though in practice m can be arbitrary but the length of x should be patched
     zeros to 2's power so that the bit reversal operation to the x,X's indices
     can function well. */
  assert((m&(m-1))==0);
  assert(m<=(fft->L));
  for (size_t i=0; i<m; i++)
    X[i] = x[mpi_br_width(i, ntdk_log2(m))];
  size_t width = ntdk_get_nbits(m);
  for (size_t logm=1; logm<width; logm++)
  {
    for (size_t j=0; j<m; j+=(1<<logm))
    {
      size_t half = 1<<(logm-1);
      for (size_t i=0; i<half; i++)
      {
        size_t idx_even = j+i;
        size_t idx_odd  = j+i+half;
        size_t zeta_idx = (i * fft->L)>>logm;
        _Complex double omega = fft->zetas[zeta_idx]*X[idx_odd];
        _Complex double butterfly_plus  = X[idx_even]+omega;
        _Complex double butterfly_minus = X[idx_even]-omega;
        X[idx_even] = butterfly_plus;
        X[idx_odd]  = butterfly_minus;
      }
    }
  }
}

void
fft_inv(size_t m, _Complex double *x, const _Complex double *X, const fft_ctx *fft)
{
  /* m is not the length of fft_ctx->L. It should satisfy m <= 2*poly_degree.
     Though in practice m can be arbitrary but the length of x should be patched
     zeros to 2's power so that the bit reversal operation to the x,X's indices
     can function well. */
  assert((m&(m-1))==0);
  assert(m<=(fft->L));
  for (size_t i=0; i<m; i++)
    x[i] = X[mpi_br_width(i, ntdk_log2(m))];
  size_t width = ntdk_get_nbits(m);
  for (size_t logm=1; logm<width; logm++)
  {
    for (size_t j=0; j<m; j+=(1<<logm))
    {
      size_t half = 1<<(logm-1);
      for (size_t i=0; i<half; i++)
      {
        size_t idx_even = j+i;
        size_t idx_odd  = j+i+half;
        size_t zeta_idx = (i * fft->L)>>logm;
        _Complex double omega = conj(fft->zetas[zeta_idx])*x[idx_odd];
        _Complex double butterfly_plus  = x[idx_even]+omega;
        _Complex double butterfly_minus = x[idx_even]-omega;
        x[idx_even] = butterfly_plus;
        x[idx_odd]  = butterfly_minus;
      }
    }
  }
  for (size_t i=0; i<m; i++)
    x[i] /= m;
}

void
emb_fwd(size_t m, _Complex double *X, const _Complex double *x, const fft_ctx *fft)
{
  /* m is the length of The message. It should satisfy m <= poly_degree/2.
     Though in practice m can be arbitrary but the length of x should be patched
     zeros to 2's power so that the bit reversal operation to the x,X's indices
     can function well. */
  assert((m&(m-1))==0);
  assert(m<=(fft->L)/4);
  for (size_t i=0; i<m; i++)
    X[i] = x[mpi_br_width(i, ntdk_log2(m))];
  size_t width = ntdk_get_nbits(m);
  for (size_t logm=1; logm<width; logm++)
  {
    size_t idx_mod = 1<<(logm+2);
    size_t gap = fft->L/idx_mod;
    for (size_t j=0; j<m; j+=(1<<logm))
    {
      size_t half = 1<<(logm-1);
      for (size_t i=0; i<half; i++)
      {
        size_t idx_even = j+i;
        size_t idx_odd  = j+i+half;
        size_t zeta_idx = (fft->rot_group[i]%idx_mod)*gap;
        _Complex double omega = fft->zetas[zeta_idx]*X[idx_odd];
        _Complex double butterfly_plus  = X[idx_even]+omega;
        _Complex double butterfly_minus = X[idx_even]-omega;
        X[idx_even] = butterfly_plus;
        X[idx_odd]  = butterfly_minus;
      }
    }
  }
}

void
emb_inv(size_t m, _Complex double *x, const _Complex double *X, const fft_ctx *fft)
{
  /* m is the length of The message. It should satisfy m <= poly_degree/2.
     Though in practice m can be arbitrary but the length of x should be patched
     zeros to 2's power so that the bit reversal operation to the x,X's indices
     can function well. */
  assert((m&(m-1))==0);
  assert(m<=(fft->L)/4);
  _Complex double *tmp = (_Complex double *)calloc(m, sizeof(_Complex double));
  for (size_t i=0; i<m; i++)
    tmp[i] = X[i];
  size_t width = ntdk_get_nbits(m);
  for (size_t logm=width-1; logm>0; logm--)
  {
    size_t idx_mod = 1<<(logm+2);
    size_t gap = fft->L/idx_mod;
    for (size_t j=0; j<m; j+=(1<<logm))
    {
      size_t half = 1<<(logm-1);
      for (size_t i=0; i<half; i++)
      {
        size_t idx_even = j+i;
        size_t idx_odd  = j+i+half;
        size_t zeta_idx = (fft->rot_group[i]%idx_mod)*gap;
        _Complex double butterfly_plus  = tmp[idx_even] + tmp[idx_odd];
        _Complex double butterfly_minus = tmp[idx_even] - tmp[idx_odd];
        butterfly_minus *= conj(fft->zetas[zeta_idx]);
        tmp[idx_even] = butterfly_plus;
        tmp[idx_odd]  = butterfly_minus;
      }
    }
  }
  for (size_t i=0; i<m; i++)
    x[i] = tmp[mpi_br_width(i, ntdk_log2(m))]/m;
}

void
poly_mod(MPI *r, const MPI Q)
{
  for (size_t i=0; i<GPQHE_N; i++)
    gcry_mpi_mod(r[i], r[i], Q);
}

void
poly_smod(MPI *r, const MPI Q)
{
  MPI Q_half = gcry_mpi_new(0); /* Q_half = Q/2 */
  gcry_mpi_div(Q_half, NULL, Q, GCRYMPI_CONST_TWO, -1);
  for (size_t i=0; i<GPQHE_N; i++)
  {
    gcry_mpi_mod(r[i], r[i], Q);          /* r[i] = r[i]%Q (same as poly_mod) */
    if (gcry_mpi_cmp(r[i], Q_half)>=0) /* if r[i] >= Q/2 */
      gcry_mpi_sub(r[i], r[i], Q);        /* r[i] = r[i]%Q -Q */
  }
  gcry_mpi_release(Q_half);
}
/*/
void
poly_round(size_t d, MPI *r, const polyd *x)
{
  for (size_t i=0; i<d; i++)
  {
    int64_t val;
    if (x->coeffs[i]>=0)
      val = (int64_t)(x->coeffs[i]+0.5);
    else
      val = (int64_t)(x->coeffs[i]-0.5);
    r->coeffs[i] = val;
  }
}
/*/
void
poly_xpy(MPI *r,
  const MPI *x, const MPI *y,
  const MPI Q, unsigned int flag)
{
  assert((flag==0)||(flag==1));
  for(size_t i=0; i<GPQHE_N; i++)
    gcry_mpi_add(r[i], x[i], y[i]);
  if (Q)
  {
    if (flag==0)
      poly_mod(r, Q);
    else if (flag==1)
      poly_smod(r, Q);
  }
}

void
poly_xsy(MPI *r,
  const MPI *x, const MPI *y,
  const MPI Q, unsigned int flag)
{
  assert((flag==0)||(flag==1));
  for(size_t i=0; i<GPQHE_N; i++)
    gcry_mpi_sub(r[i], x[i], y[i]);
  if (Q)
  {
    if (flag==0)
      poly_mod(r, Q);
    else if (flag==1)
      poly_smod(r, Q);
  }
}

void
poly_xmy(MPI *r,
  const MPI *x, const MPI *y,
  const MPI Q, unsigned int flag)
{
  assert((flag==0)||(flag==1));
  MPI *xy = gcry_calloc(GPQHE_N, sizeof(MPI));
  for (size_t k=0; k<2*GPQHE_N-1; k++) /* k = 0,...,d-1; d,...,2d-2 */
  {
    /* Since x^d = -1, the degree is taken mod d, and the sign changes when the
       exponent is > d. */
    MPI ck  = gcry_mpi_set_ui(NULL, 0);
    MPI tmp = gcry_mpi_new(0);
    for (size_t i=0; i<GPQHE_N; i++)
    {
      size_t j = k-i;
      if ((0<=j) && (j<GPQHE_N))
      {
        gcry_mpi_mul(tmp, x[i], y[j]);
        gcry_mpi_add(ck, ck, tmp);    /* ck += x[i]*y[j]; */
      }
    }
    if (k<GPQHE_N)
    {
      xy[k] = gcry_mpi_set_ui(NULL, 0);
      gcry_mpi_add(xy[k%GPQHE_N], xy[k%GPQHE_N], ck);
    }
    else
      gcry_mpi_sub(xy[k%GPQHE_N], xy[k%GPQHE_N], ck);
    gcry_mpi_release(tmp);
    gcry_mpi_release(ck);
  }
  for (size_t i=0; i<GPQHE_N; i++)
    r[i] = gcry_mpi_copy(xy[i]);
  gcry_free(xy);
  if (Q)
  {
    if (flag==0)
      poly_mod(r, Q);
    else if (flag==1)
      poly_smod(r, Q);
  }
}

void
poly_xmy_crt(size_t d, MPI *r, const MPI *x, const MPI *y, crt_ctx *crt, ntt_ctx *ntts)
{
  MPI *CRT_rep = gcry_calloc(crt->L, d*sizeof(MPI));
  MPI *values  = gcry_calloc(crt->L,   sizeof(MPI));
  for (size_t i=0; i<crt->L; i++)
  {
    values[i] = gcry_mpi_new(0);
    MPI *X  = gcry_calloc(d, sizeof(MPI));
    MPI *Y  = gcry_calloc(d, sizeof(MPI));
    MPI *XY = gcry_calloc(d, sizeof(MPI));
    for (size_t j=0; j<d; j++)
    {
      X[j]  = gcry_mpi_new(0);
      Y[j]  = gcry_mpi_new(0);
      XY[j] = gcry_mpi_new(0);
      CRT_rep[d*i+j] = gcry_mpi_new(0);
    }
    ntt_fwd(d, X, x, crt->p[i], &ntts[i]);
    ntt_fwd(d, Y, y, crt->p[i], &ntts[i]);
    for (size_t j=0; j<d; j++)
      gcry_mpi_mul(XY[j], X[j], Y[j]);
    ntt_inv(d, CRT_rep+i*d, XY, crt->p[i], &ntts[i]);
  }
  for (size_t j=0; j<d; j++)
  {
    for (size_t p=0; p<crt->L; p++)
      values[p] = gcry_mpi_copy(CRT_rep[p*d+j]);
    crt_inv(&r[j], values, crt);
    ntdk_smod(&r[j], crt->M);
  }
}

void
poly_xmy_fft(polyd *r, const polyd *x, const polyd *y)
{
/*
  assert(((r->d)==(x->d)) && ((r->d)==(y->d)));
  size_t d = r->d;
  fft_ctx fft;
  fft_build(&fft, d*2);
  _Complex double *X  = (_Complex double *)calloc(d*2, sizeof(_Complex double));
  _Complex double *Y  = (_Complex double *)calloc(d*2, sizeof(_Complex double));
  _Complex double *XY = (_Complex double *)calloc(d*2, sizeof(_Complex double));
  _Complex double *xy = (_Complex double *)calloc(d*2, sizeof(_Complex double));
  for (size_t i=0; i<d*2; i++)
  {
    X[i] = x->coeffs[i];
    Y[i] = y->coeffs[i];
  }
  fft_fwd(d, X, x->coeffs, &fft);
  fft_fwd(d, Y, y->coeffs, &fft);
  for (size_t i=0; i<fft->L; i++)
    XY[i]= X[i]*Y[i];
  int64_t *values = (int64_t *)calloc(d, sizeof(int64_t));
  for (size_t i=0; i<2*d-1; i++)
  {
    size_t idx = i%d;
    ssize_t sign = (int(i<d)-0.5)*2;
    values[idx]+=sign*
  }
*/
}

void
poly_scal(size_t d, MPI *r, const long alpha)
{
  MPI mpi_alpha = gcry_int64_to_mpi(alpha);
  for (size_t i=0;i<d; i++)
    gcry_mpi_mul(r[i], mpi_alpha, r[i]);
  gcry_mpi_release(mpi_alpha);
}

void
poly_rot(MPI *r, const MPI *x, size_t rot)
{
  size_t power = 1;
  size_t m = 2*GPQHE_N;
  for (size_t j=0; j<rot; j++)
    power *= 5; /* power = pow(5,rot) */
  for (size_t i=0; i<GPQHE_N; i++)
  {
    size_t k = (i*power)%m;
    if (k<GPQHE_N)
      r[k] = gcry_mpi_copy(x[i]); /* r[k] = x[i] */
    else
      gcry_mpi_neg(r[k-GPQHE_N], x[i]); /* r[k-d] = -x[i] */
  }
}

void
poly_conj(MPI *r, const MPI *x)
{
  r[0] = gcry_mpi_copy(x[0]);
  for (size_t i=1; i<GPQHE_N; i++)
    gcry_mpi_neg(r[i], x[GPQHE_N-i]);
}

MPI
poly_eval(size_t d, MPI *r, const int64_t x0)
{
  MPI ret = gcry_mpi_copy(r[d-1]);
//  int64_t ret = r->coeffs[r->d-1];
  for (int i=d-2; i>-1; i--)
  {
    gcry_mpi_mul(ret, ret, gcry_int64_to_mpi(x0));
    gcry_mpi_add(ret, ret, r[i]);
    //ret = ret*x0 + r->coeffs[i];
  }
  return ret;
}

/* end ntdk.c */
