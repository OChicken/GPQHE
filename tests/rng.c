/*
 * Test file: usage of rng.
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

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#define KYBER_SECRETKEYBYTES 2400
#define CRYPTO_SECRETKEYBYTES  KYBER_SECRETKEYBYTES

extern void randombytes(uint8_t *out, size_t outlen);
#ifdef AES256
extern void randombytes_init(uint8_t *entropy_input, uint8_t *personalization_string);
#endif

static uint32_t bytes_to_u32(uint8_t *x, size_t xlen)
{
  uint32_t y=0;
  for (size_t i=0; i<xlen; i++)
    y+=x[i]<<(8*i);
  return y;
}

int main()
{
#ifdef AES256 /* this is used in NIST PQC submission */
  unsigned char entropy_input[48];
  for (int i=0; i<48; i++)
    entropy_input[i] = i;
  randombytes_init(entropy_input, NULL);
#endif

  uint8_t sk[CRYPTO_SECRETKEYBYTES];
  randombytes(sk, CRYPTO_SECRETKEYBYTES);
  for (size_t i=0; i<CRYPTO_SECRETKEYBYTES; i++)
    printf("%02x", sk[i]);
  puts("");

  uint8_t bytes[sizeof(uint32_t)];
  uint32_t num;
  randombytes(bytes, sizeof(uint32_t));
  memcpy(&num, bytes, sizeof(uint32_t));
  assert(num==bytes_to_u32(bytes, sizeof(uint32_t)));

  for (size_t i=0; i<sizeof(uint32_t); i++)
    printf("%02x", bytes[i]);
  printf("=%u\n", num);

  printf("bytes[0]<< 0 = %i\n", bytes[0]<<0);
  printf("bytes[1]<< 8 = %i\n", bytes[1]<<8);
  printf("bytes[2]<<16 = %i\n", bytes[2]<<16);
  printf("bytes[3]<<24 = %i\n", bytes[3]<<24);
  printf("num=%u\n", bytes_to_u32(bytes, sizeof(uint32_t)));
  return 0;
}
