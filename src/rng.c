/*
 * rng - random number generator.
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

#include "config.h"
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h> /* abort */
#include <errno.h>
#include <fcntl.h>   /* open & creat, return int */
#include <unistd.h>  /* read & write, return ssize_t */
#include <gcrypt.h>

BEGIN_DECLS

#ifdef SUPERCOP
/* Deterministic randombytes by Daniel J. Bernstein */
/* taken from SUPERCOP (https://bench.cr.yp.to)     */

static uint32_t seed[32] = {
  3,1,4,1,5,9,2,6,5,3,5,8,9,7,9,3,2,3,8,4,6,2,6,4,3,3,8,3,2,7,9,5
};
static uint32_t in[12];
static uint32_t out[8];
static int outleft = 0;

#define ROTATE(x,b) (((x) << (b)) | ((x) >> (32 - (b))))
#define MUSH(i,b) x = t[i] += (((x ^ seed[i]) + sum) ^ ROTATE(x,b));

static void surf(void)
{
  uint32_t t[12]; uint32_t x; uint32_t sum = 0;
  int r; int i; int loop;

  for (i = 0;i < 12;++i) t[i] = in[i] ^ seed[12 + i];
  for (i = 0;i < 8;++i) out[i] = seed[24 + i];
  x = t[11];
  for (loop = 0;loop < 2;++loop) {
    for (r = 0;r < 16;++r) {
      sum += 0x9e3779b9;
      MUSH(0,5) MUSH(1,7) MUSH(2,9) MUSH(3,13)
      MUSH(4,5) MUSH(5,7) MUSH(6,9) MUSH(7,13)
      MUSH(8,5) MUSH(9,7) MUSH(10,9) MUSH(11,13)
    }
    for (i = 0;i < 8;++i) out[i] ^= t[i + 4];
  }
}

void randombytes(uint8_t *x,size_t xlen)
{
  while (xlen > 0) {
    if (!outleft) {
      if (!++in[0]) if (!++in[1]) if (!++in[2]) ++in[3];
      surf();
      outleft = 8;
    }
    *x = out[--outleft];
    ++x;
    --xlen;
  }
}
#endif

#ifdef RANDOM
void randombytes(uint8_t *out, size_t outlen)
{
  static int fd = -1;
  ssize_t ret;

  while(fd == -1) {
    fd = open("/dev/urandom", O_RDONLY);
    if(fd == -1 && errno == EINTR)
      continue;
    else if(fd == -1)
      abort();
  }

  while(outlen > 0) {
    ret = read(fd, out, outlen);
    if(ret == -1 && errno == EINTR)
      continue;
    else if(ret == -1)
      abort();

    out += ret;
    outlen -= ret;
  }
}
#endif

#ifdef GCRY_RANDOM
void randombytes(uint8_t *x, size_t xlen)
{
  gcry_randomize(x, xlen, GCRY_STRONG_RANDOM);
}
#endif

#ifdef AES256

typedef struct {
    uint8_t   Key[32];
    uint8_t   V[16];
    int       reseed_counter;
} AES256_CTR_DRBG_struct;

static AES256_CTR_DRBG_struct DRBG_ctx;

/**
 * AES256_ECB - symmetric encryption using AES256 algorithm with ECB mode
 * 
 * This uses AES from Libgcrypt.
 * 
 * @key: 256-bit AES key
 * @ctr: a 128-bit plaintext value
 * @buffer: a 128-bit ciphertext value
 */
static void AES256_ECB(uint8_t *key, uint8_t *ctr, uint8_t *buffer)
{
  gcry_cipher_hd_t hd;
  gcry_error_t err = 0;

  /* Create and initialise the context */
  err = gcry_cipher_open (&hd, GCRY_CIPHER_AES256, GCRY_CIPHER_MODE_ECB, 0);
  if (err)
    fprintf(stderr, "aes256-ecb, gcry_cipher_open failed: %s\n", gpg_strerror(err));
  err = gcry_cipher_setkey (hd, key,
    gcry_cipher_get_algo_keylen(GCRY_CIPHER_AES256));
  if (err)
    fprintf(stderr, "aes256-ecb, gcry_cipher_setkey failed: %s\n", gpg_strerror(err));
  err = gcry_cipher_encrypt (hd, buffer, 16, ctr, 16);
  if (err)
    fprintf(stderr, "aes256-ecb, gcry_cipher_encrypt failed: %s\n", gpg_strerror(err));

  /* Clean up */
  gcry_cipher_close (hd);
  if (err)
    abort();
}

static void AES256_CTR_DRBG_Update(uint8_t *provided_data, uint8_t *Key, uint8_t *V)
{
  uint8_t temp[48];

  for (int i=0; i<3; i++)
  {
    //increment V
    for (int j=15; j>=0; j--)
    {
      if ( V[j] == 0xff )
        V[j] = 0x00;
      else
      {
        V[j]++;
        break;
      }
    }
    AES256_ECB(Key, V, temp+16*i);
  }
  if ( provided_data != NULL )
    for (int i=0; i<48; i++)
      temp[i] ^= provided_data[i];
  memcpy(Key, temp, 32);
  memcpy(V, temp+32, 16);
}

void randombytes_init(uint8_t *entropy_input, uint8_t *personalization_string)
{
  uint8_t seed_material[48];

  memcpy(seed_material, entropy_input, 48);
  if (personalization_string)
    for (int i=0; i<48; i++)
      seed_material[i] ^= personalization_string[i];
  memset(DRBG_ctx.Key, 0x00, 32);
  memset(DRBG_ctx.V, 0x00, 16);
  AES256_CTR_DRBG_Update(seed_material, DRBG_ctx.Key, DRBG_ctx.V);
  DRBG_ctx.reseed_counter = 1;
}

void randombytes(uint8_t *out, size_t outlen)
{
  uint8_t  block[16];
  int      i = 0;

  while ( outlen > 0 ) {
    //increment V
    for (int j=15; j>=0; j--) {
      if ( DRBG_ctx.V[j] == 0xff )
        DRBG_ctx.V[j] = 0x00;
      else {
        DRBG_ctx.V[j]++;
        break;
      }
    }
    AES256_ECB(DRBG_ctx.Key, DRBG_ctx.V, block);
    if ( outlen > 15 ) {
      memcpy(out+i, block, 16);
      i += 16;
      outlen -= 16;
    }
    else {
      memcpy(out+i, block, outlen);
      outlen = 0;
    }
  }
  AES256_CTR_DRBG_Update(NULL, DRBG_ctx.Key, DRBG_ctx.V);
  DRBG_ctx.reseed_counter++;
}

#endif

END_DECLS
