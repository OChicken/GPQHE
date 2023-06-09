/*
 * Configurations.
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

#ifndef CONFIG_H
#define CONFIG_H

/* C++ needs to know that types and declarations are C, not C++.  */
#ifdef  __cplusplus
# define BEGIN_DECLS  extern "C" {
# define END_DECLS  }
#else
# define BEGIN_DECLS
# define END_DECLS
#endif

/* The size of `unsigned short', as computed by sizeof. */
_Static_assert(sizeof(unsigned short) == 2,
  "Require sizeof(unsigned short) == 2");
#define SIZEOF_UNSIGNED_SHORT 2

/* The size of `unsigned int', as computed by sizeof. */
_Static_assert(sizeof(unsigned int) == 4,
  "Require sizeof(unsigned int) == 4");
#define SIZEOF_UNSIGNED_INT 4

/* The size of `unsigned long', as computed by sizeof. */
_Static_assert(sizeof(unsigned long) == 8,
  "Require sizeof(unsigned long) == 8");
#define SIZEOF_UNSIGNED_LONG 8

/* The size of `unsigned long long', as computed by sizeof. */
_Static_assert(sizeof(unsigned long long) == 8,
  "Require sizeof(unsigned long long) == 8");
#define SIZEOF_UNSIGNED_LONG_LONG 8

/* The size of `void *', as computed by sizeof. */
_Static_assert(sizeof(void *) == 8,
  "Require sizeof(void *) == 8");
#define SIZEOF_VOID_P 8

#include <pthread.h>
#define GPQHE_NUM_THREAD 4

#include <assert.h>

#endif /* CONFIG_H */
