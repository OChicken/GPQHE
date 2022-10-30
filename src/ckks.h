/*
 * CKKS scheme implementation
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

/* begin ckks.h */

#ifndef CKKS_H
#define CKKS_H

#include "ntdk.h"
#include "config.h"

#define GPQHE_CKKS_SLOT GPQHE_N/2
#define GPQHE_SYMBYTES GPQHE_N//32   /* size in bytes of hashes, and seeds */
#define GPQHE_INDCPA_MSGBYTES       GPQHE_SYMBYTES
#define GPQHE_INDCPA_SECRETKEYBYTES GPQHE_N
#define GPQHE_INDCPA_BYTES GPQHE_N
#define RotGroup 5

BEGIN_DECLS

struct ckks_ctx{
  size_t slots;
  double Delta;
  MPI ql_table[GPQHE_L];
  double Bclean;                   /* Freshly encrypted noise upper bound */
  double Brs;                      /* rescale noise upper bound */
  double Bks;                      /* key switch noise upper bound */
  crt_ctx crt;
  size_t pbnd;                     /* `pbnd=59` in HEAAN */
  size_t num_p;                    /* the number of primes for CRT context */
  ntt_ctx *ntts;
  _Complex double *CRT_U;
  _Complex double *CRT_UConj;
  _Complex double *CRT_UTrans;
  _Complex double *CRT_UConjTrans;
  double cM;                       /* ring constant depending on M(=2d) */
};

struct ckks_ct{
  size_t l;                        /* Current evaluation level, 1<=l<=L */
  double nu;
  MPI c0[GPQHE_INDCPA_BYTES];
  MPI c1[GPQHE_INDCPA_BYTES];
};

struct ckks_pk{
  MPI p0[GPQHE_INDCPA_BYTES];
  MPI p1[GPQHE_INDCPA_BYTES];
};

/** Build CKKS context */
void he_ckks_ctx_build(int precomp);

/** Show CKKS context parameters */
void he_ckks_ctx_show_params(void);

/** Show CKKS ciphertext parameters */
void he_ckks_ct_show_params(const struct ckks_ct *ct, const char *title, ...);

/*******************************************************************************
 * Encode                                                                      *
 ******************************************************************************/

/** Encodes complex array into pt, an integral polynomial. */
void he_ckks_encode(double pt[GPQHE_INDCPA_MSGBYTES],
                    const size_t m,
                    const COMPLEX *msg);

/** Decodes pt, an integral polynomial, into complex array. */
void he_ckks_decode(size_t n, COMPLEX *msg, double pt[GPQHE_INDCPA_MSGBYTES]);

/** Encodes a constant to pt */
void he_ckks_const_pt(double pt[GPQHE_INDCPA_MSGBYTES], COMPLEX num);

/** Canonical embedding norm */
double he_ckks_canemb_norm(const double pt[GPQHE_INDCPA_MSGBYTES]);

/*******************************************************************************
 * Encryption/Decryption                                                       *
 ******************************************************************************/

/** Generate sk and pk */
void indcpa_keypair(int8_t sk[GPQHE_INDCPA_SECRETKEYBYTES],
                    struct ckks_pk *pk);

/** Generate relinearization key */
void he_ckks_genrlk(struct ckks_pk *rlk,
                    const int8_t sk[GPQHE_N]);

/** Generate conjugate key */
void he_ckks_genck(struct ckks_pk *ck,
                   const int8_t sk[GPQHE_N]);

/** Generate rotation key */
void he_ckks_genrk(struct ckks_pk *rk,
                   const int8_t sk[GPQHE_N]);

/** public key encryption */
void indcpa_enc_pk(struct ckks_ct *ct,
                   const double pt[GPQHE_INDCPA_MSGBYTES],
                   const struct ckks_pk *pk);
/*
void indcpa_enc(uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                const uint8_t coins[KYBER_SYMBYTES])
*/

/** secret key encryption */
void indcpa_enc_sk(struct ckks_ct *ct,
                   const double pt[GPQHE_INDCPA_MSGBYTES],
                   const int8_t sk[GPQHE_INDCPA_SECRETKEYBYTES]);

/** decryption */
void indcpa_dec(double pt[GPQHE_INDCPA_MSGBYTES],
                const struct ckks_ct *ct,
                const int8_t sk[GPQHE_INDCPA_SECRETKEYBYTES]);
/*
void indcpa_dec(uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES])
*/

void he_ckks_ct_copy(struct ckks_ct *ct_dest, const struct ckks_ct *ct);

/*******************************************************************************
 * Homomorphic Operations                                                      *
 ******************************************************************************/

/** ct + ct */
void he_ckks_add(struct ckks_ct *ct_dest,
                 const struct ckks_ct *ct1,
                 const struct ckks_ct *ct2);

/** ct + pt */
void he_ckks_add_pt(struct ckks_ct *ct_dest,
                    const struct ckks_ct *ct,
                    const double pt[GPQHE_INDCPA_MSGBYTES]);

/** -ct */
void he_ckks_neg(struct ckks_ct *ct);

/** ct - ct */
void he_ckks_sub(struct ckks_ct *ct_dest,
  const struct ckks_ct *ct1, const struct ckks_ct *ct2);

/** ct - pt */
void he_ckks_sub_pt(struct ckks_ct *ct_dest,
                    const struct ckks_ct *ct,
                    const double pt[GPQHE_INDCPA_MSGBYTES]);

/** ct * ct */
void he_ckks_mult(struct ckks_ct *ct_dest,
                  const struct ckks_ct *ct1,
                  const struct ckks_ct *ct2,
                  const struct ckks_pk *rlk);

/** ct * pt */
void he_ckks_mult_pt(struct ckks_ct *ct_dest,
                     const struct ckks_ct *ct,
                     const double pt[GPQHE_INDCPA_MSGBYTES]);

/** conjugate */
void he_ckks_conj(struct ckks_ct *ct_dest,
                  const struct ckks_ct *ct,
                  const struct ckks_pk *ck);

/** rotation */
void he_ckks_rot(struct ckks_ct *ct_dest,
                 const struct ckks_ct *ct,
                 const size_t rot,
                 const struct ckks_pk *rk);

/** rescale */
void he_ckks_rescale(struct ckks_ct *ct);

/** mod down */
void he_ckks_moddown(struct ckks_ct *ct);

/*******************************************************************************
 * Advanced                                                                    *
 ******************************************************************************/

/** gemv: GEneral Matrix Vector multiplication
 * he_ckks_gemv(&ct_Av, A, &ct_v, evk->rot_key, ctx); */
void he_ckks_gemv(struct ckks_ct *ct_dest,
                  double _Complex *mat,
                  const struct ckks_ct *ct,
                  const struct ckks_pk *rk);

/** coefficient to slots
 * Takes an encryption of t(x) = t_0 + t_1x + ... and transforms to encryptions
 * of (t_0, t_1, ..., t_(d/2)) and (t_(d/2 + 1), ..., t_(n-1)) before these
 * vectors are encoded. */
void he_ckks_coeff_to_slot(struct ckks_ct *ct_real,
                           struct ckks_ct *ct_imag,
                           const struct ckks_ct *ct,
                           const struct ckks_pk *ck,
                           const struct ckks_pk *rk);

/** slots to coefficients
 * Takes encryptions of (t_0, t_1, ..., t_(n/2-1)) and (t_(n/2), ..., t_(n-1))
 * before these vectors are encoded and transofmrs to an encryption of
 * t(x) = t_0 + t_1x + ... */
void he_ckks_slot_to_coeff(struct ckks_ct *ct,
                           const struct ckks_ct *ct0,
                           const struct ckks_ct *ct1,
                           const struct ckks_pk *rk);

/** Exaluates the exponential function on the ciphertext. Takes an encryption
 * of ct and a and returns exp(a*ct). */
void he_ckks_exp(struct ckks_ct *ct_dest,
                 const struct ckks_ct *ct,
                 const struct ckks_pk *rlk);

/** Exaluates the sinusoidal function on the ciphertext. Takes an encryption
 * of ct and a and returns sin(a*ct). */
void he_ckks_sin(struct ckks_ct *ct_dest,
                 const struct ckks_ct *ct,
                 const struct ckks_pk *rlk);

void he_ckks_rlsin(struct ckks_ct *ct_dest,
                   const struct ckks_ct *ct,
                   const double a,
                   const struct ckks_pk *rlk,
                   const struct ckks_pk *ck);

void he_ckks_bootstrap(struct ckks_ct *ct,
                       const struct ckks_pk *rlk,
                       const struct ckks_pk *ck,
                       const struct ckks_pk *rk);

/** inv(ct) - do element-wise inverse 1/ct */
void he_ckks_inv(struct ckks_ct *ct_inv,
                 const struct ckks_ct *ct,
                 const struct ckks_pk *rlk,
                 const unsigned int iter);

/** sqrt(ct) */
void he_ckks_sqrt(struct ckks_ct *ct_sqrt,
                  const struct ckks_ct *ct,
                  const struct ckks_pk *rlk,
                  const unsigned int iter);

END_DECLS

#endif /* CKKS_H */
