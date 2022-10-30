/* 
 * CKKS scheme implementation
 * Copyright (C) shouran.ma@rwth-aachen.de
 *
 * This file is part of GPQHE.
 * GPQHE is designed on top of NewHope, Kyber and HEAAN.
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

#include "ckks.h"

#include <fcntl.h>   /* open & creat, return int */
#include <unistd.h>  /* read & write, return ssize_t */

#include <pmu.h>

#define UINT64MAX_PER_DBLMAX 15

#ifndef GPQHE_PRECOMP
#define GPQHE_PRECOMP 1
#endif

//static struct ring_ctx ring;
static struct ckks_ctx ctx;

/* Build CKKS context */
void he_ckks_ctx_build(int precomp)
{
  /* initialize */
#if 0
  if (precomp)
  {
    ring_ctx_build(&ring);
    ring_ctx_save(&ring);
  }
  else
    ring_ctx_load(&ring);
#endif
  gcry_control (GCRYCTL_DISABLE_SECMEM, 0);
  gcry_control (GCRYCTL_ENABLE_QUICK_RANDOM, 0);
  gcry_control (GCRYCTL_INITIALIZATION_FINISHED, 0);

  size_t d = GPQHE_N;//ring.n;
  ctx.slots = d/2;

  /* Scaling factor */
  ctx.Delta = GPQHE_P;

  /* modulus ql = q0*p^l */
  for (size_t l=0; l<GPQHE_L; l++)
  {
    ctx.ql_table[l] = gcry_mpi_new(0);
    if (l==0)
      mpi_set_ui(ctx.ql_table[0], GPQHE_P);//gcry_mpi_set_ui(ctx.ql_table[0], ring.p_vec[0]);
    else
#if 0
      gcry_mpi_mul_ui(ctx.ql_table[l], ctx.ql_table[l-1], ring.p_vec[l]);
#else
      gcry_mpi_mul_ui(ctx.ql_table[l], ctx.ql_table[l-1], GPQHE_P);
#endif
  }

  /* bounds */
  size_t h = GPQHE_HW; /* hamming weight */
  double sigma = GPQHE_SIGMA;
  ctx.Bclean = 8*sqrt(2)*sigma*d + 6*sigma*sqrt(d) + 16*sigma*sqrt(h*d);
  ctx.Brs = sqrt((double)d/3) * (3+8*sqrt(h));
  ctx.Bks = 8*sigma*d/sqrt(3);
  assert(GPQHE_P > (GPQHE_N + 2*ctx.Bclean));

#if 0
  /* crt build */
  ctx.pbnd = 59;
  MPI qL = he_modulus_ql(ctx.L, ctx);
  unsigned int log2_qL = gcry_mpi_get_nbits(qL)-1;
  ctx.num_p = (2 + ntdk_log2(GPQHE_N) + 4*log2_qL + ctx.pbnd - 1) / ctx.pbnd;
  crt_build(&ctx.crt, ctx.pbnd, ctx.num_p);
#endif

  /* ntt build */
  //ctx.ntts = (ntt_ctx *)calloc(ctx.pbnd, sizeof(ntt_ctx));
  //for (size_t i=0; i<ctx.num_p; i++)
    //ntt_build(&ctx.ntts[i], ctx.d, ctx.crt.p[i], ctx.crt.g[i]);

  /* CRT matrix */
  size_t d_half = d/2;
  size_t      m = d*2;
  size_t  power = 1;
  size_t rot    = RotGroup;
  ctx.CRT_U          = (_Complex double *)calloc(d_half*d, sizeof(_Complex double));
  ctx.CRT_UConj      = (_Complex double *)calloc(d_half*d, sizeof(_Complex double));
  ctx.CRT_UTrans     = (_Complex double *)calloc(d*d_half, sizeof(_Complex double));
  ctx.CRT_UConjTrans = (_Complex double *)calloc(d*d_half, sizeof(_Complex double));
  double *cM_list     = (double *)calloc(d, sizeof(double));
  for (size_t i=0; i<d_half; i++)
  {
    double theta = 2*M_PI*power / m;
    _Complex double zeta = cos(theta) + I*sin(theta);
    power = (power*rot)%m;
    for (size_t j=0; j<d; j++)
    {
      if (!j)
        ctx.CRT_U[i*d+j] = 1;
      else
        ctx.CRT_U[i*d+j] = ctx.CRT_U[i*d+j-1]*zeta;
      ctx.CRT_UConj[i*d+j]           = conj(ctx.CRT_U[i*d+j]);
      ctx.CRT_UConjTrans[j*d_half+i] = conj(ctx.CRT_U[i*d+j]);
      ctx.CRT_UTrans[j*d_half+i]     =      ctx.CRT_U[i*d+j] ;
      cM_list[j] += (cabs(ctx.CRT_UConjTrans[j*d_half+i])
                 +   cabs(ctx.CRT_UTrans[j*d_half+i]))/d;
    }
  }
  ctx.cM = cM_list[0];
  for (size_t i=1; i<d; i++)
  {
    if (ctx.cM < cM_list[i])
      ctx.cM = cM_list[i];
  }

  /* release */
  //gcry_mpi_release(qL);
  free(cM_list);
}

void he_ckks_ctx_show_params()
{
  printf("\033[1mEncryption parameters:\033[0m\n");
#if GPQHE_CQ == 'C'
  printf("\\_ Security type: classical %u bits.\n", GPQHE_SEC_LEVEL);
#elif GPQHE_CQ == 'Q'
  printf("\\_ Security type: quantum %u bits.\n", GPQHE_SEC_LEVEL);
#endif
  printf("\\_ Polynomial degree d = %u\n", GPQHE_N);
  printf("\\_ Cyclotomic index M = %u\n", 2*GPQHE_N);
  printf("\\_ Ring constant cM = %f\n", ctx.cM);
  printf("\\_ Level L = %i\n", GPQHE_L);

  /* Q = q1*q1*...*qL */
  printf("\\_ Ciphertext modulus ql = q0*p^l,\n");
  printf("   \\_ q0 (%u bits) = ", gcry_mpi_get_nbits(ctx.ql_table[0]));
  show_mpi(ctx.ql_table[0]);
  printf("   \\_ p = 1<<%lu = %#02llX\n", ntdk_log2(GPQHE_P), GPQHE_P);
  printf("   \\_ qL (%u bits) = ", gcry_mpi_get_nbits(ctx.ql_table[GPQHE_L-1]));
  show_mpi(ctx.ql_table[GPQHE_L-1]);
  printf("\\_ Bclean = %f\n", ctx.Bclean);
  printf("\\_ Brs    = %f\n", ctx.Brs);
  printf("\\_ Bks    = %f\n", ctx.Bks);
  //printf("\\_ Number of primes for CRT context: %lu\n", ctx->num_p);
}

void he_ckks_ct_show_params(const struct ckks_ct *ct, const char *title, ...)
{
  va_list arg_ptr;
  printf("\033[1m\033[4m");
  va_start (arg_ptr, title) ;
  vfprintf (stdout, title, arg_ptr);
  va_end (arg_ptr);
  printf("\033[0m\n");
  printf("\\_ Level l = %lu\n", ct->l);
  printf("\\_ nu = %.2f (%ld bits)\n", ct->nu, (long)log2(ct->nu)+1);
}

/*******************************************************************************
 * Encode                                                                      *
 ******************************************************************************/

/** Encodes complex array into pt, an integral polynomial. */
void he_ckks_encode(double pt[GPQHE_INDCPA_MSGBYTES],
                    const size_t m,
                    const _Complex double *msg)
{
  assert(m<=GPQHE_N/2);
  COMPLEX msg_round[GPQHE_N/2] = {0};
  COMPLEX tmp[GPQHE_N] = {0};
  for (size_t i=0; i<m; i++)
  {
    int sign_real = creal(msg[i])>=0? 1 : -1;
    int sign_imag = cimag(msg[i])>=0? 1 : -1;
    msg_round[i]  =   (double)((int64_t)(creal(msg[i])*ctx.Delta + sign_real*0.5))
                  + I*(double)((int64_t)(cimag(msg[i])*ctx.Delta + sign_imag*0.5));
  }
  for (size_t i=0; i<GPQHE_N; i++)
  {
    tmp[i] = 0;
    for (size_t j=0; j<GPQHE_N/2; j++)
    {
      if (j<m)
        tmp[i] += (ctx.CRT_UConjTrans[i*GPQHE_N/2+j]*     msg_round[j]
               +   ctx.CRT_UTrans    [i*GPQHE_N/2+j]*conj(msg_round[j]))/GPQHE_N;
      else
        tmp[i] += 0;
    }
    assert(cimag(tmp[i]/ctx.Delta)<FLT_EPSILON);
    pt[i] = creal(tmp[i]);
  }
}

/** Decodes pt, an integral polynomial, into complex array. */
void
he_ckks_decode(size_t n, _Complex double *msg, double pt[GPQHE_INDCPA_MSGBYTES])
{
  assert(n<=(GPQHE_N/2));
  for (size_t i=0; i<n; i++)
  {
    msg[i] = 0;
    for (size_t j=0;j<GPQHE_N;j++)
      msg[i] += ctx.CRT_U[i*GPQHE_N+j]*pt[j]/ctx.Delta;
  }
}

/** Encodes a constant to pt */
void he_ckks_const_pt(double pt[GPQHE_INDCPA_MSGBYTES], COMPLEX num)
{
  for (size_t i=0; i<GPQHE_N; i++)
    pt[i]=0;
  pt[0]=num*GPQHE_P;
}

/** Canonical embedding norm */
double he_ckks_canemb_norm(const double pt[GPQHE_INDCPA_MSGBYTES])
{
  size_t d_half = GPQHE_N/2;
  double canemb_norm = 0;
  for (size_t i=0; i<d_half; i++)
  {
    _Complex double val=0;
    for (size_t j=0; j<GPQHE_N; j++)
      val += ctx.CRT_U[i*GPQHE_N+j]*pt[j];
    if (!i)
      canemb_norm = cabs(val);
    else
    {
      if (canemb_norm < cabs(val))
        canemb_norm = cabs(val);
    }
  }
  return canemb_norm;
}

/*******************************************************************************
 * Encryption/Decryption                                                       *
 ******************************************************************************/

/** Generate keys */
static void
he_genkey_switching_key(struct ckks_pk *swk,
  const MPI *new_key, const MPI *sk,
  const MPI qL, const MPI P)
{
  MPI *e   = gcry_calloc(GPQHE_N, sizeof(MPI));
  int16_t e_int16[GPQHE_N] = {0};
  sample_discrete_gaussian(GPQHE_N, e_int16, GPQHE_SIGMA);
  MPI PqL = gcry_mpi_new(0);
  gcry_mpi_mul(PqL, P, qL);
  unsigned int nbits_PqL = gcry_mpi_get_nbits(PqL);
  for (size_t i=0; i<GPQHE_N; i++)
    gcry_mpi_mul(new_key[i], new_key[i], P);
  for (size_t i=0; i<GPQHE_N; i++)
  {
    e[i] = gcry_int64_to_mpi((int64_t)e_int16[i]);
    do {
      gcry_mpi_randomize(swk->p1[i], nbits_PqL, GCRY_WEAK_RANDOM);
    } while(gcry_mpi_cmp(swk->p1[i], PqL)>=0);
    /* p1[i] is as long as PqL but less than P.qL */
  }
  poly_xmy(swk->p0, swk->p1, sk, NULL, 0);
  for (size_t i=0; i<GPQHE_N; i++)
    gcry_mpi_neg(swk->p0[i], swk->p0[i]);
  poly_xpy(swk->p0, swk->p0, e, NULL, 0);      /* -sw1.sk + e */
  poly_xpy(swk->p0, swk->p0, new_key, PqL, 1); /* (-sw1.sk + e + P.newkey) mod PqL */
  gcry_free(e);
  gcry_mpi_release(PqL);
}

void indcpa_keypair(int8_t sk[GPQHE_INDCPA_SECRETKEYBYTES],
                    struct ckks_pk *pk)
{
  printf("Generating sk and pk ... ");
  fflush(stdout);
  /* secret key */
  sample_hwt(GPQHE_N, sk, GPQHE_HW);
  /* sample error for public key */
  int16_t e_int16[GPQHE_N];
  sample_discrete_gaussian(GPQHE_N, e_int16, GPQHE_SIGMA);
  /* pk initialization */
  MPI ql = gcry_mpi_copy(ctx.ql_table[GPQHE_L-1]);
  MPI e[GPQHE_N];
  MPI s[GPQHE_N];
  unsigned int nbits_qL = gcry_mpi_get_nbits(ql);
  for (size_t i=0; i<GPQHE_N; i++){
    e[i]      = gcry_int64_to_mpi((int64_t)e_int16[i]);
    s[i]      = gcry_int64_to_mpi((int64_t)sk[i]);
    pk->p0[i] = gcry_mpi_new(0);
    pk->p1[i] = gcry_mpi_new(0);
    gcry_mpi_randomize(pk->p1[i], nbits_qL, GCRY_WEAK_RANDOM);
    /* p1[i] is as long as qL but less than qL */
    while(gcry_mpi_cmp(pk->p1[i], ql)>=0)
      gcry_mpi_sub(pk->p1[i], pk->p1[i], ql);
  }
  poly_xmy(pk->p0, pk->p1, s, NULL, 0);
  for (size_t i=0; i<GPQHE_N; i++)
    gcry_mpi_neg(pk->p0[i], pk->p0[i]);
  poly_xpy(pk->p0, pk->p0, e, ql, 1);
  /* release */
  gcry_mpi_release(ql);
  for (size_t i=0; i<GPQHE_N; i++)
  {
    gcry_mpi_release(e[i]);
    gcry_mpi_release(s[i]);
  }
  printf("done.\n\n");
}

void he_ckks_genrlk(struct ckks_pk *rlk,
                    const int8_t sk[GPQHE_N])
{
  printf("Generating rlk ... ");
  fflush(stdout);
  MPI ql = gcry_mpi_copy(ctx.ql_table[GPQHE_L-1]);
  MPI P  = gcry_mpi_copy(ql);
  MPI s [GPQHE_N];
  MPI s2[GPQHE_N];
  for (size_t i=0; i<GPQHE_N; i++)
  {
    s [i]      = gcry_int64_to_mpi((int64_t)sk[i]);
    s2[i]      = gcry_mpi_new(0);
    rlk->p0[i] = gcry_mpi_new(0);
    rlk->p1[i] = gcry_mpi_new(0);
  }
  /* relinearization key */
  poly_xmy(s2, s, s, NULL, 0);
  he_genkey_switching_key(rlk, s2, s, ql, P);
  /* release */
  gcry_mpi_release(ql);
  gcry_mpi_release(P);
  for (size_t i=0; i<GPQHE_N; i++)
  {
    gcry_mpi_release(s [i]);
    gcry_mpi_release(s2[i]);
  }
  printf("done.\n\n");
}

void he_ckks_genck(struct ckks_pk *ck,
                   const int8_t sk[GPQHE_N])
{
  printf("Generating ck ... ");
  fflush(stdout);
  MPI ql = gcry_mpi_copy(ctx.ql_table[GPQHE_L-1]);
  MPI P  = gcry_mpi_copy(ql);
  MPI  s [GPQHE_N];
  MPI _ck[GPQHE_N];
  for (size_t i=0; i<GPQHE_N; i++)
  {
     s [i]      = gcry_int64_to_mpi((int64_t)sk[i]);
    _ck[i]      = gcry_mpi_new(0);
    ck->p0[i] = gcry_mpi_new(0);
    ck->p1[i] = gcry_mpi_new(0);
  }
  /* conjugate key */
  poly_conj(_ck, s);
  he_genkey_switching_key(ck, _ck, s, ql, P);
  /* release */
  gcry_mpi_release(ql);
  gcry_mpi_release(P);
  for (size_t i=0; i<GPQHE_N; i++)
  {
    gcry_mpi_release( s [i]);
    gcry_mpi_release(_ck[i]);
  }
  printf("done.\n\n");
}

void he_ckks_genrk(struct ckks_pk *rk,
                   const int8_t sk[GPQHE_N])
{
  printf("Generating rk ... ");
  fflush(stdout);
  MPI ql = gcry_mpi_copy(ctx.ql_table[GPQHE_L-1]);
  MPI P  = gcry_mpi_copy(ql);
  MPI  s [GPQHE_N];
  MPI _rk[GPQHE_N];
  for (size_t i=0; i<GPQHE_N; i++)
  {
     s [i]     = gcry_int64_to_mpi((int64_t)sk[i]);
    _rk[i]     = gcry_mpi_new(0);
    for (size_t rot=0; rot<GPQHE_CKKS_SLOT; rot++)
    {
      rk[rot].p0[i] = gcry_mpi_new(0);
      rk[rot].p1[i] = gcry_mpi_new(0);
    }
  }
  /* rotation keys */
  for (size_t rot=0; rot<GPQHE_CKKS_SLOT; rot++)
  {
    poly_rot(_rk, s, rot);
    he_genkey_switching_key(&rk[rot], _rk, s, ql, P);
  }
  /* release */
  gcry_mpi_release(ql);
  gcry_mpi_release(P);
  for (size_t i=0; i<GPQHE_N; i++)
  {
    gcry_mpi_release( s [i]);
    gcry_mpi_release(_rk[i]);
  }
  printf("done.\n\n");
}

void indcpa_enc_pk(struct ckks_ct *ct,
                   const double pt[GPQHE_INDCPA_MSGBYTES],
                   const struct ckks_pk *pk)
{
  ct->l = GPQHE_L;
  ct->nu = ctx.Delta;
  MPI ql = gcry_mpi_copy(ctx.ql_table[ct->l-1]);
  int16_t v_int16[GPQHE_N];
  int16_t e0_int16[GPQHE_N];
  int16_t e1_int16[GPQHE_N];
  sample_zero_center(GPQHE_N, v_int16);
  sample_discrete_gaussian(GPQHE_N, e0_int16, GPQHE_SIGMA);
  sample_discrete_gaussian(GPQHE_N, e1_int16, GPQHE_SIGMA);
  MPI m [GPQHE_N];
  MPI v [GPQHE_N];
  MPI e0[GPQHE_N];
  MPI e1[GPQHE_N];
  for (size_t i=0; i<GPQHE_N; i++)
  {
    m[i]      = gcry_int128_to_mpi((__int128_t)pt[i]);
    v[i]      = gcry_int64_to_mpi((int64_t)v_int16[i]);
    e0[i]     = gcry_int64_to_mpi((int64_t)e0_int16[i]);
    e1[i]     = gcry_int64_to_mpi((int64_t)e1_int16[i]);
    ct->c0[i] = gcry_mpi_new(0);
    ct->c1[i] = gcry_mpi_new(0);
  }
  poly_xmy(ct->c0, pk->p0, v , NULL, 0);
  poly_xpy(ct->c0, ct->c0, m , NULL, 0);
  poly_xpy(ct->c0, ct->c0, e0,   ql, 1); /* (v*pk + pt + e0) smod qL */
  poly_xmy(ct->c1, pk->p1, v , NULL, 0);
  poly_xpy(ct->c1, ct->c1, e1,   ql, 1); /* (v*pk      + e1) smod qL */
  /* release */
  gcry_mpi_release(ql);
  for (size_t i=0; i<GPQHE_N; i++)
  {
    gcry_mpi_release(m [i]);
    gcry_mpi_release(v [i]);
    gcry_mpi_release(e0[i]);
    gcry_mpi_release(e1[i]);
  }
}

void indcpa_enc_sk(struct ckks_ct *ct,
                   const double pt[GPQHE_INDCPA_MSGBYTES],
                   const int8_t sk[GPQHE_INDCPA_SECRETKEYBYTES])
{
  ct->l = GPQHE_L;
  ct->nu = ctx.Delta;
  MPI ql = gcry_mpi_copy(ctx.ql_table[ct->l-1]);
  unsigned int nbits_qL = gcry_mpi_get_nbits(ql);
  int16_t v_int16[GPQHE_N];
  int16_t e_int16[GPQHE_N];
  sample_zero_center(GPQHE_N, v_int16);
  sample_discrete_gaussian(GPQHE_N, e_int16, GPQHE_SIGMA);
  MPI m[GPQHE_N];
  MPI v[GPQHE_N];
  MPI s[GPQHE_N];
  MPI e[GPQHE_N];
  for (size_t i=0; i<GPQHE_N; i++)
  {
    m[i]      = gcry_int128_to_mpi((__int128_t)pt[i]);
    v[i]      = gcry_int64_to_mpi((int64_t)v_int16[i]);
    s[i]      = gcry_int64_to_mpi((int64_t)sk[i]);
    e[i]      = gcry_int64_to_mpi((int64_t)e_int16[i]);
    ct->c0[i] = gcry_mpi_new(0);
    ct->c1[i] = gcry_mpi_new(0);
    gcry_mpi_randomize(ct->c1[i], nbits_qL, GCRY_WEAK_RANDOM);
    while(gcry_mpi_cmp(ct->c1[i], ql)>=0)
      gcry_mpi_sub(ct->c1[i], ct->c1[i], ql);
  }
  poly_xmy(ct->c0, ct->c1, s, NULL, 0);
  for (size_t i=0; i<GPQHE_N; i++)
    gcry_mpi_neg(ct->c0[i], ct->c0[i]);
  poly_xpy(ct->c0, ct->c0, m, NULL, 0);
  poly_xpy(ct->c0, ct->c0, e,   ql, 1); /* (-c1*sk + pt + e) smod qL */
  poly_smod(ct->c1, ql);                /* c1 smod qL */
  /* release */
  gcry_mpi_release(ql);
  for (size_t i=0; i<GPQHE_N; i++)
  {
    gcry_mpi_release(m[i]);
    gcry_mpi_release(v[i]);
    gcry_mpi_release(s[i]);
    gcry_mpi_release(e[i]);
  }
}

void indcpa_dec(double pt[GPQHE_INDCPA_MSGBYTES],
                const struct ckks_ct *ct,
                const int8_t sk[GPQHE_INDCPA_SECRETKEYBYTES])
{
  MPI ql = gcry_mpi_copy(ctx.ql_table[ct->l-1]);
  MPI m[GPQHE_N];
  MPI s[GPQHE_N];
  for (size_t i=0; i<GPQHE_N; i++)
  {
    m[i] = gcry_mpi_new(0);
    s[i] = gcry_int64_to_mpi((int64_t)sk[i]);
  }
  poly_xmy(m, ct->c1, s, NULL, 0);
  poly_xpy(m, ct->c0, m,   ql, 1); /* (c0 + c1*sk) smod ql */
  for (size_t i=0; i<GPQHE_N; i++)
    pt[i] = (double)gcry_mpi_to_int128(m[i]);
  /* release */
  gcry_mpi_release(ql);
  for (size_t i=0; i<GPQHE_N; i++)
  {
    gcry_mpi_release(m[i]);
    gcry_mpi_release(s[i]);
  }
}

/*******************************************************************************
 * Copy Plaintext & Ciphertext (internal used)                                 *
 ******************************************************************************/

/** copy ciphertext */
void
he_ckks_ct_copy(struct ckks_ct *ct_dest, const struct ckks_ct *ct)
{
  ct_dest->l  = ct->l;
  ct_dest->nu = ct->nu;
  for (size_t i=0; i<GPQHE_N; i++)
  {
    ct_dest->c0[i] = gcry_mpi_copy(ct->c0[i]);
    ct_dest->c1[i] = gcry_mpi_copy(ct->c1[i]);
  }
}

/*******************************************************************************
 * Refresh (internal used)                                                     *
 ******************************************************************************/

/** Relinearizes a dim3 ciphertext back down to dim2 */
static void
he_ckks_relinearize(struct ckks_ct *ct_dest,
  const MPI *d0, const MPI *d1, const MPI *d2,
  const struct ckks_pk *rlk, const MPI ql, const MPI P)
{
  /* "We note that P·qL is the largest modulus to generate evaluation key and it
   * suffices to assume that P is approximately equal to qL." --- CKKS 2017 */
  poly_xmy(ct_dest->c0, d2, rlk->p0, NULL, 0);
  poly_xmy(ct_dest->c1, d2, rlk->p1, NULL, 0);
  for (size_t i=0; i<GPQHE_N; i++)
  {
    ntdk_div_round(ct_dest->c0[i], ct_dest->c0[i], P);
    ntdk_div_round(ct_dest->c1[i], ct_dest->c1[i], P);
  }
  poly_xpy(ct_dest->c0, ct_dest->c0, d0, ql, 1);
  poly_xpy(ct_dest->c1, ct_dest->c1, d1, ql, 1);
}

/** Key-switching procedure */
static void
he_ckks_switch_key(struct ckks_ct *ct_dest,
  const MPI *d0, const MPI *d1,
  const struct ckks_pk *swk, const MPI ql, const MPI P)
{
  MPI Pql = gcry_mpi_new(0);
  gcry_mpi_mul(Pql, P, ql);
  poly_xmy(ct_dest->c0, d1, swk->p0, Pql, 1);
  poly_xmy(ct_dest->c1, d1, swk->p1, Pql, 1);
  for (size_t i=0; i<GPQHE_N; i++)
  {
    ntdk_div_round(ct_dest->c0[i], ct_dest->c0[i], P);
    ntdk_div_round(ct_dest->c1[i], ct_dest->c1[i], P);
  }
  poly_xpy (ct_dest->c0, ct_dest->c0, d0, ql, 1);
  poly_smod(ct_dest->c1, ql);
  gcry_mpi_release(Pql);
}

/*******************************************************************************
 * Homomorphic Operations                                                      *
 ******************************************************************************/

/** ct_dest = ct1+ct2 */
void
he_ckks_add(struct ckks_ct *ct_dest,
            const struct ckks_ct *ct1,
            const struct ckks_ct *ct2)
{
  assert(ct1->l==ct2->l);
  if ((ct_dest!=ct1)&&(ct_dest!=ct2))
    he_ckks_ct_copy(ct_dest, ct1);
  ct_dest->nu = ct1->nu;//+ct2->nu;
  MPI ql = gcry_mpi_copy(ctx.ql_table[ct_dest->l-1]);
  poly_xpy(ct_dest->c0, ct1->c0, ct2->c0, ql, 1);
  poly_xpy(ct_dest->c1, ct1->c1, ct2->c1, ql, 1);
  gcry_mpi_release(ql);
}

/** ct_dest = ct+pt */
void
he_ckks_add_pt(struct ckks_ct *ct_dest,
               const struct ckks_ct *ct,
               const double pt[GPQHE_INDCPA_MSGBYTES])
{
  if (ct_dest!=ct)
    he_ckks_ct_copy(ct_dest, ct);
  ct_dest->l  = ct->l;
  ct_dest->nu = ct->nu;
  MPI ql = gcry_mpi_copy(ctx.ql_table[ct_dest->l-1]);
  MPI m[GPQHE_N];
  for (size_t i=0; i<GPQHE_N; i++)
    m[i] = gcry_int128_to_mpi((__int128_t)pt[i]);
  poly_xpy (ct_dest->c0, ct->c0, m, ql, 1);
  poly_smod(ct_dest->c1, ql);
  /* release */
  gcry_mpi_release(ql);
  for (size_t i=0; i<GPQHE_N; i++)
    gcry_mpi_release(m[i]);
}

/** -ct */
void
he_ckks_neg(struct ckks_ct *ct)
{
  MPI ql = gcry_mpi_copy(ctx.ql_table[ct->l-1]);
  MPI ql_half = gcry_mpi_new(0);
  gcry_mpi_div(ql_half, NULL, ql, GCRYMPI_CONST_TWO, -1);
  for (size_t i=0; i<GPQHE_N; i++)
  {
    gcry_mpi_neg(ct->c0[i], ct->c0[i]);
    gcry_mpi_neg(ct->c1[i], ct->c1[i]);
    gcry_mpi_mod(ct->c0[i], ct->c0[i], ql);
    gcry_mpi_mod(ct->c1[i], ct->c1[i], ql);
    if (gcry_mpi_cmp(ct->c0[i], ql_half)>=0)
      gcry_mpi_sub(ct->c0[i], ct->c0[i], ql);
    if (gcry_mpi_cmp(ct->c1[i], ql_half)>=0)
      gcry_mpi_sub(ct->c1[i], ct->c1[i], ql);
  }
  /* release */
  gcry_mpi_release(ql);
  gcry_mpi_release(ql_half);
}

/** ct_dest = ct1-ct2 */
void
he_ckks_sub(struct ckks_ct *ct_dest,
  const struct ckks_ct *ct1, const struct ckks_ct *ct2)
{
  assert(ct1->l==ct2->l);
  if ((ct_dest!=ct1)&&(ct_dest!=ct2))
    he_ckks_ct_copy(ct_dest, ct1);
  MPI ql = gcry_mpi_copy(ctx.ql_table[ct_dest->l-1]);
  ct_dest->nu = ct1->nu;//+ct2->nu;
  poly_xsy(ct_dest->c0, ct1->c0, ct2->c0, ql, 1);
  poly_xsy(ct_dest->c1, ct1->c1, ct2->c1, ql, 1);
  /* release */
  gcry_mpi_release(ql);
}

/** ct_dest = ct-pt */
void
he_ckks_sub_pt(struct ckks_ct *ct_dest,
               const struct ckks_ct *ct,
               const double pt[GPQHE_INDCPA_MSGBYTES])
{
  if (ct_dest!=ct)
    he_ckks_ct_copy(ct_dest, ct);
  ct_dest->l  = ct->l;
  ct_dest->nu = ct->nu;
  MPI ql = gcry_mpi_copy(ctx.ql_table[ct_dest->l-1]);
  MPI m[GPQHE_N];
  for (size_t i=0; i<GPQHE_N; i++)
    m[i] = gcry_int128_to_mpi((__int128_t)pt[i]);
  poly_xsy(ct_dest->c0, ct->c0, m, ql, 1);
  poly_smod(ct_dest->c1, ql);
  /* release */
  gcry_mpi_release(ql);
  for (size_t i=0; i<GPQHE_N; i++)
    gcry_mpi_release(m[i]);
}

/** ct_dest = ct1*ct2 */
void he_ckks_mult(struct ckks_ct *ct_dest,
                  const struct ckks_ct *ct1,
                  const struct ckks_ct *ct2,
                  const struct ckks_pk *rlk)
{
  assert(ct1->l==ct2->l);
  if ((ct_dest!=ct1)&&(ct_dest!=ct2))
    he_ckks_ct_copy(ct_dest, ct1);
  ct_dest->nu = ct1->nu*ct2->nu;
  MPI ql = gcry_mpi_copy(ctx.ql_table[ct_dest->l-1]);
  MPI P    = gcry_mpi_copy(ql);
  MPI d0 [GPQHE_N];
  MPI d1 [GPQHE_N];
  MPI d11[GPQHE_N];
  MPI d12[GPQHE_N];
  MPI d2 [GPQHE_N];
  for (size_t i=0; i<GPQHE_N; i++)
  {
    d0[i]  = gcry_mpi_new(0);
    d1[i]  = gcry_mpi_new(0);
    d11[i] = gcry_mpi_new(0);
    d12[i] = gcry_mpi_new(0);
    d2[i]  = gcry_mpi_new(0);
  }
  poly_xmy(d0, ct1->c0, ct2->c0,    ql, 1);
  poly_xmy(d2, ct1->c1, ct2->c1,    ql, 1);
  poly_xmy(d11, ct1->c0, ct2->c1, NULL, 0);
  poly_xmy(d12, ct1->c1, ct2->c0, NULL, 0);
  poly_xpy(d1, d11, d12, ql, 1);
  he_ckks_relinearize(ct_dest, d0, d1, d2, rlk, ql, P);
  /* release */
  gcry_mpi_release(ql);
  gcry_mpi_release(P);
  for (size_t i=0; i<GPQHE_N; i++)
  {
    gcry_mpi_release(d0 [i]);
    gcry_mpi_release(d1 [i]);
    gcry_mpi_release(d11[i]);
    gcry_mpi_release(d12[i]);
    gcry_mpi_release(d2 [i]);
  }
}

/** ct_dest = ct*pt */
void he_ckks_mult_pt(struct ckks_ct *ct_dest,
                     const struct ckks_ct *ct,
                     const double pt[GPQHE_INDCPA_MSGBYTES])
{
  if (ct_dest!=ct)
    he_ckks_ct_copy(ct_dest, ct);
  ct_dest->nu = ct->nu*ct->nu; /* FIX ME */
  MPI ql = gcry_mpi_copy(ctx.ql_table[ct_dest->l-1]);
  MPI m[GPQHE_N];
  for (size_t i=0; i<GPQHE_N; i++)
    m[i] = gcry_int128_to_mpi((__int128_t)pt[i]);
  poly_xmy(ct_dest->c0, m, ct->c0, ql, 1);
  poly_xmy(ct_dest->c1, m, ct->c1, ql, 1);
  /* release */
  gcry_mpi_release(ql);
  for (size_t i=0; i<GPQHE_N; i++)
    gcry_mpi_release(m[i]);
}

void he_ckks_conj(struct ckks_ct *ct_dest,
                  const struct ckks_ct *ct,
                  const struct ckks_pk *ck)
{
  he_ckks_ct_copy(ct_dest, ct);
  MPI ql = gcry_mpi_copy(ctx.ql_table[ct_dest->l-1]);
  MPI P  = gcry_mpi_copy(ql);
  MPI d0[GPQHE_N];
  MPI d1[GPQHE_N];
  for (size_t i=0; i<GPQHE_N; i++)
  {
    d0[i] = gcry_mpi_new(0);
    d1[i] = gcry_mpi_new(0);
  }
  poly_conj(d0, ct->c0);
  poly_conj(d1, ct->c1);
  he_ckks_switch_key(ct_dest, d0, d1, ck, ql, P);
  /* release */
  gcry_mpi_release(ql);
  gcry_mpi_release(P);
  for (size_t i=0; i<GPQHE_N; i++)
  {
    gcry_mpi_release(d0[i]);
    gcry_mpi_release(d1[i]);
  }
}

void he_ckks_rot(struct ckks_ct *ct_dest,
                 const struct ckks_ct *ct,
                 const size_t rot,
                 const struct ckks_pk *rk)
{
  assert(ct_dest!=ct);
  he_ckks_ct_copy(ct_dest, ct);
  MPI ql = gcry_mpi_copy(ctx.ql_table[ct_dest->l-1]);
  MPI P  = gcry_mpi_copy(ql);
  MPI d0[GPQHE_N];
  MPI d1[GPQHE_N];
  for (size_t i=0; i<GPQHE_N; i++)
  {
    d0[i] = gcry_mpi_new(0);
    d1[i] = gcry_mpi_new(0);
  }
  poly_rot(d0, ct->c0, rot);
  poly_rot(d1, ct->c1, rot);
  he_ckks_switch_key(ct_dest, d0, d1, rk, ql, P);
  /* release */
  gcry_mpi_release(ql);
  gcry_mpi_release(P);
  for (size_t i=0; i<GPQHE_N; i++)
  {
    gcry_mpi_release(d0[i]);
    gcry_mpi_release(d1[i]);
  }
}

/** rescale */
void
he_ckks_rescale(struct ckks_ct *ct)
{
#if 0
  MPI p = gcry_int64_to_mpi(ring.p_vec[ct->l-1]);
#else
  MPI p = gcry_int64_to_mpi(GPQHE_P);
#endif
  ct->l -= 1;
  ct->nu /= ctx.Delta;
  MPI ql = gcry_mpi_copy(ctx.ql_table[ct->l-1]);
  MPI ql_half = gcry_mpi_new(0);
  gcry_mpi_div(ql_half, NULL, ql, GCRYMPI_CONST_TWO, -1);
  for (size_t i=0; i<GPQHE_N; i++)
  {
    ntdk_div_round(ct->c0[i],ct->c0[i],p);
    ntdk_div_round(ct->c1[i],ct->c1[i],p);
    /* signed mod */
    gcry_mpi_mod(ct->c0[i], ct->c0[i], ql);
    if (gcry_mpi_cmp(ct->c0[i], ql_half)>=0)
      gcry_mpi_sub(ct->c0[i], ct->c0[i], ql);
    gcry_mpi_mod(ct->c1[i], ct->c1[i], ql);
    if (gcry_mpi_cmp(ct->c1[i], ql_half)>=0)
      gcry_mpi_sub(ct->c1[i], ct->c1[i], ql);
  }
  gcry_mpi_release(p);
  gcry_mpi_release(ql);
  gcry_mpi_release(ql_half);
}

/** mod down */
void he_ckks_moddown(struct ckks_ct *ct)
{
  ct->l -= 1;
  MPI ql = gcry_mpi_copy(ctx.ql_table[ct->l-1]);
  MPI ql_half = gcry_mpi_new(0);
  gcry_mpi_div(ql_half, NULL, ql, GCRYMPI_CONST_TWO, -1);
  for (size_t i=0; i<GPQHE_N; i++)
  {
    gcry_mpi_mod(ct->c0[i], ct->c0[i], ql);
    if (gcry_mpi_cmp(ct->c0[i], ql_half)>=0)
      gcry_mpi_sub(ct->c0[i], ct->c0[i], ql);
    gcry_mpi_mod(ct->c1[i], ct->c1[i], ql);
    if (gcry_mpi_cmp(ct->c1[i], ql_half)>=0)
      gcry_mpi_sub(ct->c1[i], ct->c1[i], ql);
  }
  gcry_mpi_release(ql);
  gcry_mpi_release(ql_half);
}

/*******************************************************************************
 * Advanced                                                                    *
 ******************************************************************************/

/** gemv: GEneral Matrix Vector multiplication
 * he_ckks_gemv(&ct_Av, A, &ct_v, evk->rot_key, ctx); */
void he_ckks_gemv(struct ckks_ct *ct_dest,
                  double _Complex *mat,
                  const struct ckks_ct *ct,
                  const struct ckks_pk *rk)
{
  /* plaintext mat is of dim (d/2)*(d/2).  */
  size_t d_half = GPQHE_N/2;
  size_t N1 = ntdk_sqrt_u64(d_half);
  if (d_half != N1*N1)
    N1 = ntdk_sqrt_u64(2*d_half);
  size_t N2 = d_half/N1;

  /* Compute sum */
  struct ckks_ct outer_sum;
  for (size_t i=0; i<N2; i++)
  {
    size_t shift = i*N1;
    struct ckks_ct inner_sum;
    for (size_t j=0; j<N1; j++)
    {
      /* rotate and encode ct */
      struct ckks_ct ct_rot;
      he_ckks_rot(&ct_rot, ct, j, &rk[j]);
      /* rotate and encode mat */
      COMPLEX diag    [GPQHE_N/2] = {0};
      COMPLEX diag_rot[GPQHE_N/2] = {0};
      double pt[GPQHE_INDCPA_MSGBYTES];
      blas_zdiag(d_half, diag, shift+j, mat);
      blas_zrot(d_half, diag_rot, (int64_t)(-shift), diag);
      he_ckks_encode(pt, d_half, diag_rot);
      /* multiply and accumulate */
      struct ckks_ct dot_prod;
      he_ckks_mult_pt(&dot_prod, &ct_rot, pt);
      if (j)
        he_ckks_add(&inner_sum, &inner_sum, &dot_prod);
      else
        he_ckks_ct_copy(&inner_sum, &dot_prod);
    }
    struct ckks_ct rotated_sum;
    he_ckks_rot(&rotated_sum, &inner_sum, shift, &rk[shift]);
    if (i)
      he_ckks_add(&outer_sum, &outer_sum, &rotated_sum);
    else
      he_ckks_ct_copy(&outer_sum, &rotated_sum);
  }
  he_ckks_ct_copy(ct_dest, &outer_sum);
  he_ckks_rescale(ct_dest);
}


/** coefficient to slots
 * Takes an encryption of t(x) = t_0 + t_1x + ... and transforms to encryptions
 * of (t_0, t_1, ..., t_(d/2)) and (t_(d/2 + 1), ..., t_(n-1)) before these
 * vectors are encoded. */
void he_ckks_coeff_to_slot(struct ckks_ct *ct_real,
                           struct ckks_ct *ct_imag,
                           const struct ckks_ct *ct,
                           const struct ckks_pk *ck,
                           const struct ckks_pk *rk)
{
  COMPLEX U0_ConjTrans[GPQHE_CKKS_SLOT*GPQHE_CKKS_SLOT];
  COMPLEX U0_Trans[GPQHE_CKKS_SLOT*GPQHE_CKKS_SLOT];
  COMPLEX U1_ConjTrans[GPQHE_CKKS_SLOT*GPQHE_CKKS_SLOT];
  COMPLEX U1_Trans[GPQHE_CKKS_SLOT*GPQHE_CKKS_SLOT];
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
  {
    for (size_t j=0; j<GPQHE_CKKS_SLOT; j++)
    {
      U0_ConjTrans[i*GPQHE_CKKS_SLOT+j] = ctx.CRT_UConjTrans[i*GPQHE_CKKS_SLOT+j];
      U1_ConjTrans[i*GPQHE_CKKS_SLOT+j] = ctx.CRT_UConjTrans[(i+GPQHE_CKKS_SLOT)*GPQHE_CKKS_SLOT+j];
      U0_Trans[i*GPQHE_CKKS_SLOT+j] = ctx.CRT_UTrans[i*GPQHE_CKKS_SLOT+j];
      U1_Trans[i*GPQHE_CKKS_SLOT+j] = ctx.CRT_UTrans[(i+GPQHE_CKKS_SLOT)*GPQHE_CKKS_SLOT+j];
    }
  }

  struct ckks_ct ct_conj,ct0,ct1;
  double pt[GPQHE_INDCPA_MSGBYTES];
  he_ckks_const_pt(pt, (double)1/GPQHE_N);
  he_ckks_conj(&ct_conj, ct, ck);

  /* ct_real */
  he_ckks_gemv(&ct0, U0_ConjTrans,  ct     , rk);
  he_ckks_gemv(&ct1, U0_Trans    , &ct_conj, rk);
  he_ckks_add(ct_real, &ct0, &ct1);
  he_ckks_mult_pt(ct_real, ct_real, pt);
  he_ckks_rescale(ct_real);

  /* ct_imag */
  he_ckks_gemv(&ct0, U1_ConjTrans,  ct     , rk);
  he_ckks_gemv(&ct1, U1_Trans    , &ct_conj, rk);
  he_ckks_add(ct_imag, &ct0, &ct1);
  he_ckks_mult_pt(ct_imag, ct_imag, pt);
  he_ckks_rescale(ct_imag);
}

/** slots to coefficients
 * Takes encryptions of (t_0, t_1, ..., t_(n/2-1)) and (t_(n/2), ..., t_(n-1))
 * before these vectors are encoded and transofmrs to an encryption of
 * t(x) = t_0 + t_1x + ... */
void he_ckks_slot_to_coeff(struct ckks_ct *ct,
                           const struct ckks_ct *ct0,
                           const struct ckks_ct *ct1,
                           const struct ckks_pk *rk)
{
  COMPLEX U0[GPQHE_CKKS_SLOT*GPQHE_CKKS_SLOT];
  COMPLEX U1[GPQHE_CKKS_SLOT*GPQHE_CKKS_SLOT];
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
  {
    for (size_t j=0; j<GPQHE_CKKS_SLOT; j++)
    {
      U0[i*GPQHE_CKKS_SLOT+j] = ctx.CRT_U[i*GPQHE_N+j];
      U1[i*GPQHE_CKKS_SLOT+j] = ctx.CRT_U[i*GPQHE_N+j+GPQHE_CKKS_SLOT];
    }
  }
  struct ckks_ct Uct0, Uct1;
  he_ckks_gemv(&Uct0, U0, ct0, rk);
  he_ckks_gemv(&Uct1, U1, ct1, rk);
  he_ckks_add(ct, &Uct0, &Uct1);
}

/** exp(ct) */
void he_ckks_exp(struct ckks_ct *ct_dest,
                 const struct ckks_ct *ct,
                 const struct ckks_pk *rlk)
{
  double pt[GPQHE_INDCPA_MSGBYTES];
  struct ckks_ct ct2, ct4, ct01, ct23, ct0123, ct45, ct67, ct4567;
  /* ct2 = ct^2 */
  he_ckks_mult(&ct2, ct, ct, rlk);
  he_ckks_rescale(&ct2); /* l-1 */
  /* ct4 = ct^4 */
  he_ckks_mult(&ct4, &ct2, &ct2, rlk);
  he_ckks_rescale(&ct4); /* l-2 */
  /* ct01 */
  he_ckks_const_pt(pt, 1.0);
  he_ckks_add_pt(&ct01, ct, pt);
  he_ckks_mult_pt(&ct01, &ct01, pt);
  he_ckks_rescale(&ct01); /* l-1 */
  he_ckks_moddown(&ct01); /* l-2 */
  /* ct23 */
  he_ckks_const_pt(pt, 3.0);
  he_ckks_add_pt(&ct23, ct, pt);
  he_ckks_const_pt(pt, 1.0/6); /* 1/3! */
  he_ckks_mult_pt(&ct23, &ct23, pt);
  he_ckks_rescale(&ct23); /* l-1 */
  he_ckks_mult(&ct23, &ct2, &ct23, rlk);
  he_ckks_rescale(&ct23); /* l-2 */
  /* ct0123 */
  he_ckks_add(&ct0123, &ct01, &ct23); /* l-2 */
  he_ckks_moddown(&ct0123); /* l-3 */
  /* ct45 */
  he_ckks_const_pt(pt, 5.0);
  he_ckks_add_pt(&ct45, ct, pt);
  he_ckks_const_pt(pt, 1.0/120); /* 1/5! */
  he_ckks_mult_pt(&ct45, &ct45, pt);
  he_ckks_rescale(&ct45); /* l-1 */
  he_ckks_moddown(&ct45); /* l-2 */
  /* ct67 */
  he_ckks_const_pt(pt, 7.0);
  he_ckks_add_pt(&ct67, ct, pt);
  he_ckks_const_pt(pt, 1.0/5040); /* 1/7! */
  he_ckks_mult_pt(&ct67, &ct67, pt);
  he_ckks_rescale(&ct67); /* l-1 */
  he_ckks_mult(&ct67, &ct2, &ct67, rlk);
  he_ckks_rescale(&ct67); /* l-2 */
  /* ct4567 */
  he_ckks_add(&ct4567, &ct45, &ct67); /* l-2 */
  he_ckks_mult(&ct4567, &ct4, &ct4567, rlk);
  he_ckks_rescale(&ct4567); /* l-3 */
  /* ct_dest */
  he_ckks_add(ct_dest, &ct0123, &ct4567); /* l-3 */
}

#if 0
void he_ckks_sin(ckks_ct *ct_dest, const ckks_ct *ct,
  const ckks_pk *rlk, const ckks_ctx *ctx)
{
  ckks_pt pt;
  ckks_ct coeff_ct, ct2, ct4, ct13, ct57;
  /* ct2 = ct^2 */
  he_ckks_mult(&ct2, ct, ct, rlk, ctx);
  he_ckks_rescale(&ct2, ctx); /* l-1 */
  /* ct4 = ct^4 */
  he_ckks_mult(&ct4, &ct2, &ct2, rlk, ctx);
  he_ckks_rescale(&ct4, ctx); /* l-2 */
  /* ct/3! */
  he_ckks_const_pt(pt, 1.0/6);
  he_ckks_mult_pt(&coeff_ct, ct, &pt, ctx);
  he_ckks_rescale(&coeff_ct, ctx); /* x/3! (l-1) */
  /* ct13 */
  he_ckks_const_pt(&pt, 6.0);
  he_ckks_sub_pt(&ct13, &ct2, &pt, ctx); /* x^2-3! */
  he_ckks_neg(&ct13, ctx); /* 3!-x^2 (l-1) */
  he_ckks_mult(&ct13, &coeff_ct, &ct13, rlk, ctx);
  he_ckks_rescale(&ct13, ctx); /* (x/3!)(3!-x^2)=x-x^3/3! (l-2) */
  /* ct/7! */
  he_ckks_const_pt(&pt, 1.0/5040);
  he_ckks_mult_pt(&coeff_ct, ct, &pt, ctx);
  he_ckks_rescale(&coeff_ct, ctx); /* x/7! (l-1) */
  /* ct57 */
  he_ckks_const_pt(&pt, 42); /* 7!/5! = 6*7 = 42 */
  he_ckks_sub_pt(&ct57, &ct2, &pt, ctx); /* x^2-7!/5! */
  he_ckks_neg(&ct57, ctx); /* 7!/5!-x^2 (l-1) */
  he_ckks_mult(&ct57, &coeff_ct, &ct57, rlk, ctx);
  he_ckks_rescale(&ct57, ctx); /* (x/7!)(7!/5!-x^2) = x/5!-x^3/7! (l-2) */
  he_ckks_mult(&ct57, &ct57, &ct4, rlk, ctx);
  he_ckks_rescale(&ct57, ctx); /* (x/7!)(7!/5!-x^2)x^4 = x^5/5!-x^7/7! (l-3) */
  /* ct_dest */
  he_ckks_moddown(&ct13, ctx);
  he_ckks_add(ct_dest, &ct13, &ct57, ctx);
}
#endif

/**
 * he_ckks_rlsin - The relineared version of ckks_sin
 * 
 * Calculate sin(a*ct)/a = (exp(a*i*ct)-exp(-a*i*ct))/(2*i*a) 
 */ 
void he_ckks_rlsin(struct ckks_ct *ct_dest,
                   const struct ckks_ct *ct,
                   const double a,
                   const struct ckks_pk *rlk,
                   const struct ckks_pk *ck)
{
  COMPLEX const_pt[GPQHE_CKKS_SLOT];
  COMPLEX b;
  double pt[GPQHE_INDCPA_MSGBYTES];
  struct ckks_ct a_ct, ct_exp, ct_exp_neg;
  /* i*a/(2^r) */
  b = a*I/(1<<GPQHE_CKKS_ITER);
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    const_pt[i] = b;
  he_ckks_encode(pt, GPQHE_CKKS_SLOT, const_pt);
  /* i*a/(2^r)*ct */
  he_ckks_mult_pt(&a_ct, ct, pt);
  he_ckks_rescale(&a_ct); /* l-1 */
  /* exp(i*a/(2^r)*ct)^(2^r) = exp(a*i*ct) */
  he_ckks_exp(&ct_exp, &a_ct, rlk);
  for (size_t i=0; i<GPQHE_CKKS_ITER; i++)
  {
    he_ckks_mult(&ct_exp, &ct_exp, &ct_exp, rlk);
    he_ckks_rescale(&ct_exp);
  }
  /* exp(-a*i*ct) */
  he_ckks_conj(&ct_exp_neg, &ct_exp, ck);
  /* exp(a*i*ct)-exp(-a*i*ct) == 2*i*sin(a*ct) */
  he_ckks_sub(ct_dest, &ct_exp, &ct_exp_neg);
  /* 1/(2*i*a) */
  b = 1/(2*I*a);
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    const_pt[i] = b;
  he_ckks_encode(pt, GPQHE_CKKS_SLOT, const_pt);
  /* sin(a*ct)/a */
  he_ckks_mult_pt(ct_dest, ct_dest, pt);
  he_ckks_rescale(ct_dest);
}

void he_ckks_bootstrap(struct ckks_ct *ct,
                       const struct ckks_pk *rlk,
                       const struct ckks_pk *ck,
                       const struct ckks_pk *rk)
{
  /**
   * First of all, based on the well-known fact that UINT64_MAX=1.8446744e+19,
   * while DBL_MAX=1.7976931e+308, so a cascade product of 15 u64 numbers can
   * be represented by double, whereas 16 of them would cause `inf`.
   * 
   * Since the current modulus ql is the product of l limbs which do not excess
   * machine word size (u64), so l must smaller than 16, otherwise ql would be
   * represented as `inf` in double.
   */
  assert(ct->l<=UINT64MAX_PER_DBLMAX);

  double ql;
  /* Raise scaling factor. */
  double old_Delta = ctx.Delta;
  for (size_t l=0; l<ct->l; l++)
  {
    if (l==0)
      ql = GPQHE_P;//(double)ring.p_vec[0];
    else
#if 0
      ctx.Delta *= (double)ring.p_vec[l];
#else
      ql *= (double)GPQHE_P;
#endif
  }
  ctx.Delta = ql;
  ct->nu = ctx.Delta;

  /* Raises ciphertext modulus from (mod ql) to (mod qL). */
  size_t old_l = ct->l;
  ct->l = GPQHE_L;

  struct ckks_ct ct0, ct1;

  /* Coefficient to Slot */
  he_ckks_ct_show_params(ct, "before coeff_to_slot");
  he_ckks_coeff_to_slot(&ct0, &ct1, ct, ck, rk);
  he_ckks_ct_show_params(&ct0, "ct0 after coeff_to_slot, before rlsin");
  he_ckks_ct_show_params(&ct1, "ct1 after coeff_to_slot, before rlsin");

  /* Exponentiate */
  double a = 2*M_PI/ql;
  he_ckks_rlsin(&ct0, &ct0, a, rlk, ck);
  he_ckks_rlsin(&ct1, &ct1, a, rlk, ck);
  he_ckks_ct_show_params(&ct0, "ct0 after rlsin, before slot_to_coeff");
  he_ckks_ct_show_params(&ct1, "ct1 after rlsin, before slot_to_coeff");

  /* Slot to Coefficient */
  he_ckks_slot_to_coeff(ct, &ct0, &ct1, rk);
  he_ckks_ct_show_params(ct, "after slot_to_coeff");
  he_ckks_rescale(ct);
  he_ckks_ct_show_params(ct, "finally rescale");

  /* Reset scaling factor */
  ctx.Delta = old_Delta;
  ct->nu = ctx.Delta;
}

/** 1/ct
 * do inverse 1/ct element-wise */
void he_ckks_inv(struct ckks_ct *ct_inv,
                 const struct ckks_ct *ct,
                 const struct ckks_pk *rlk,
                 const unsigned int iter)
{
  double one[GPQHE_INDCPA_MSGBYTES];
  double two[GPQHE_INDCPA_MSGBYTES];
  struct ckks_ct tmp, an, bn;
  he_ckks_const_pt(one, 1);
  he_ckks_const_pt(two, 2);
  he_ckks_ct_copy(&tmp, ct);
  he_ckks_neg(&tmp);               /* tmp = -ct */
  he_ckks_add_pt(&an, &tmp, two); /* an = 2-ct */
  he_ckks_moddown(&an);
  he_ckks_add_pt(&bn, &tmp, one); /* bn = 1-ct */
  for (size_t n=0; n<iter; n++)
  {
    he_ckks_mult(&bn, &bn, &bn, rlk);  /* bn = bn*bn */
    he_ckks_rescale(&bn);
    he_ckks_add_pt(&tmp, &bn, one);         /* tmp = 1+bn */
    he_ckks_mult(&an, &an, &tmp, rlk); /* an = an*(1+bn) */
    he_ckks_rescale(&an);
  }
  he_ckks_ct_copy(ct_inv, &an);
}

void he_ckks_sqrt(struct ckks_ct *ct_sqrt,
                  const struct ckks_ct *ct,
                  const struct ckks_pk *rlk,
                  const unsigned int iter)
{
  COMPLEX num   [GPQHE_CKKS_SLOT];
  double one    [GPQHE_INDCPA_MSGBYTES];
  double three  [GPQHE_INDCPA_MSGBYTES];
  double half   [GPQHE_INDCPA_MSGBYTES];
  double quarter[GPQHE_INDCPA_MSGBYTES];
  struct ckks_ct tmp, an, bn;
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    num[i]=1;
  he_ckks_encode(one, GPQHE_CKKS_SLOT, num);
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    num[i]=3;
  he_ckks_encode(three, GPQHE_CKKS_SLOT, num);
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    num[i]=0.5;
  he_ckks_encode(half, GPQHE_CKKS_SLOT, num);
  for (size_t i=0; i<GPQHE_CKKS_SLOT; i++)
    num[i]=0.25;
  he_ckks_encode(quarter, GPQHE_CKKS_SLOT, num);
#if 0
  he_ckks_const_pt(one, 1);
  he_ckks_const_pt(three, 3);
  he_ckks_const_pt(half, 0.5);
  he_ckks_const_pt(quarter, 0.25);
#endif
  he_ckks_ct_copy(&an, ct);           /* an = ct */
  he_ckks_sub_pt(&bn, ct, one); /* bn = ct-1 */
  for (size_t n=0; n<1; n++)
  {
    /* an */
    he_ckks_mult_pt(&tmp, &bn, half); /* bn/2 */
    he_ckks_rescale(&tmp);
    he_ckks_sub_pt(&tmp, &tmp, one);  /* bn/2-1 */
    he_ckks_neg(&tmp);                 /* tmp = 1-bn/2 */
    he_ckks_moddown(&an);
    he_ckks_mult(&an, &an, &tmp, rlk); /* an = an*(1-bn/2) */
    he_ckks_rescale(&an);
    /* bn */
    he_ckks_sub_pt(&tmp, &bn, three);
    he_ckks_mult_pt(&tmp, &tmp, quarter);   /* tmp = (bn-3)/4 */
    he_ckks_rescale(&tmp);
    he_ckks_mult(&bn, &bn, &bn, rlk);  /* bn = bn*bn */
    he_ckks_rescale(&bn);
    he_ckks_mult(&bn, &bn, &tmp, rlk); /* bn = (bn*bn)*((bn-3)/4) */
    he_ckks_rescale(&bn);
  }
  he_ckks_ct_copy(ct_sqrt, &an);
  //he_ckks_ct_copy(ct_sqrt, &bn);
}