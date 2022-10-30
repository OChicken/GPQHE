/*
 * Controller parameters precomputation.
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

#ifndef CONTROLLER_H
#define CONTROLLER_H

#include "params.h"

BEGIN_DECLS

struct controller_ctx{
  double x[GPQHE_DIMX];
  double u[GPQHE_DIMU];
  double A[GPQHE_DIMX*GPQHE_DIMX];
  double B[GPQHE_DIMX*GPQHE_DIMU];
  double Q[GPQHE_DIMX*GPQHE_DIMX];
  double R[GPQHE_DIMU*GPQHE_DIMU];
  double P[GPQHE_DIMX*GPQHE_DIMX];
  double K[(GPQHE_DIMU*GPQHE_HORIZON)*GPQHE_DIMX];
};

void linear_controller_ctx_build(struct controller_ctx *controller);
void linear_controller_ctx_save (struct controller_ctx *controller);
void linear_controller_ctx_load (struct controller_ctx *controller);

END_DECLS

#endif /* CONTROLLER_H */
