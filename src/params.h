/*
 * Parameters definitions.
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

#ifndef PARAMS_H
#define PARAMS_H

#include "types.h"

#ifndef GPQHE_LOGP
  #define GPQHE_LOGP 59u  /* each prime must be larger than 59 bits */
#endif

#ifdef TEST_CRT
  #undef GPQHE_LOGP
  #define GPQHE_LOGP 9u
#endif

/**
 * GPQHE_CQ - classical or quantum
 * Available values are 'C' (for classical) or 'Q' (for quantum)
 */
#define GPQHE_CQ 'C'

/*
 * GPQHE_SEC_LEVEL - security level
 * Change this for different security strengths.
 * Available values are 128, 192, 256
 */
#define GPQHE_SEC_LEVEL 128u

#define GPQHE_K 8
#define GPQHE_BLKSIZ 64   /* block size */
#define GPQHE_SYMBYTES 32 /* size in bytes of hashes, and seeds */

#define GPQHE_PI 3.141592653589793238462643383279502884

#define GPQHE_ROT 5
#define GPQHE_SIGMA (double)3.1915382432114616

#endif /* PARAMS_H */
