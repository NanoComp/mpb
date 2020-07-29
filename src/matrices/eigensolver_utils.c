/* Copyright (C) 1999-2020 Massachusetts Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "config.h"
#include <mpiglue.h>
#include <check.h>
#include <scalar.h>
#include <matrices.h>

#include "eigensolver.h"

/**************************************************************************/

void eigensolver_get_eigenvals_aux(evectmatrix Y, real *eigenvals, evectoperator A, void *Adata,
                                   evectmatrix Work1, evectmatrix Work2, sqmatrix U, sqmatrix Usqrt,
                                   sqmatrix Uwork) {
  sqmatrix_sqrt(Usqrt, U, Uwork);       /* Usqrt = 1/sqrt(Yt*Y) */
  evectmatrix_XeYS(Work1, Y, Usqrt, 1); /* Work1 = orthonormalize(Y) */

  A(Work1, Work2, Adata, 1, Y);            /* Work2 = A Work1; Y is scratch */
  evectmatrix_XtY(U, Work1, Work2, Uwork); /* U = Work1 * A * Work1 */

  sqmatrix_eigensolve(U, eigenvals, Uwork);
  evectmatrix_XeYS(Y, Work1, U, 1);
}

void eigensolver_get_eigenvals(evectmatrix Y, real *eigenvals, evectoperator A, void *Adata,
                               evectmatrix Work1, evectmatrix Work2) {
  sqmatrix U, Usqrt, Uwork;

  U = create_sqmatrix(Y.p);
  Usqrt = create_sqmatrix(Y.p);
  Uwork = create_sqmatrix(Y.p);

  evectmatrix_XtX(U, Y, Uwork);
  sqmatrix_invert(U, 1, Uwork);

  eigensolver_get_eigenvals_aux(Y, eigenvals, A, Adata, Work1, Work2, U, Usqrt, Uwork);

  destroy_sqmatrix(U);
  destroy_sqmatrix(Usqrt);
  destroy_sqmatrix(Uwork);
}

/**************************************************************************/

/* Subroutines for chaining constraints, to make it easy to pass
   multiple constraint functions to the eigensolver: */

evectconstraint_chain *evect_add_constraint(evectconstraint_chain *constraints, evectconstraint C,
                                            void *constraint_data) {
  evectconstraint_chain *new_constraints;

  CHK_MALLOC(new_constraints, evectconstraint_chain, 1);

  new_constraints->C = C;
  new_constraints->constraint_data = constraint_data;
  new_constraints->next = constraints;
  return new_constraints;
}

void evect_destroy_constraints(evectconstraint_chain *constraints) {
  while (constraints) {
    evectconstraint_chain *cur_constraint = constraints;
    constraints = constraints->next;
    free(cur_constraint);
  }
}

void evectconstraint_chain_func(evectmatrix X, void *data) {
  evectconstraint_chain *constraints = (evectconstraint_chain *)data;

  while (constraints) {
    if (constraints->C) constraints->C(X, constraints->constraint_data);
    constraints = constraints->next;
  }
}
