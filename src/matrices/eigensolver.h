/* Copyright (C) 1999, 2000, 2001, 2002, Massachusetts Institute of Technology.
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

#ifndef EIGENSOLVER_H
#define EIGENSOLVER_H

#include "matrices.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

typedef void (*evectoperator) (evectmatrix Xin, evectmatrix Xout,
			       void *data, int is_current_eigenvector,
			       evectmatrix Work);

typedef void (*evectpreconditioner) (evectmatrix Xin, evectmatrix Xout,
				     void *data,
				     evectmatrix Y, real *eigenvals,
				     sqmatrix YtY);

typedef void (*evectconstraint) (evectmatrix X, void *data);

extern void eigensolver(evectmatrix Y, real *eigenvals,
			evectoperator A, void *Adata,
			evectpreconditioner K, void *Kdata,
			evectconstraint constraint, void *constraint_data,
			evectmatrix Work[], int nWork,
			real tolerance, int *num_iterations,
			int flags);

extern void eigensolver_lagrange(evectmatrix Y, real *eigenvals,
			evectoperator A, void *Adata,
			evectpreconditioner K, void *Kdata,
			evectconstraint constraint, void *constraint_data,
			evectoperator L, void *Ldata, real *lag,
			evectmatrix Work[], int nWork,
			real tolerance, int *num_iterations,
			int flags);

extern void eigensolver_davidson(evectmatrix Y, real *eigenvals,
				 evectoperator A, void *Adata,
				 evectpreconditioner K, void *Kdata,
				 evectconstraint constraint,
				 void *constraint_data,
				 evectmatrix Work[], int nWork,
				 real tolerance, int *num_iterations,
				 int flags,
				 real target);

extern void eigensolver_get_eigenvals(evectmatrix Y, real *eigenvals,
				      evectoperator A, void *Adata,
				      evectmatrix Work1, evectmatrix Work2);

/* eigensolver option flags, designed to be combined with a bitwise or ('|');
   each flag should set exactly one bit. */
#define EIGS_VERBOSE (1<<0)
#define EIGS_PROJECT_PRECONDITIONING (1<<1)
#define EIGS_RESET_CG (1<<2)
#define EIGS_FORCE_EXACT_LINMIN (1<<3)
#define EIGS_FORCE_APPROX_LINMIN (1<<4)
#define EIGS_ORTHONORMALIZE_FIRST_STEP (1<<5)
#define EIGS_REORTHOGONALIZE (1<<6)
#define EIGS_DYNAMIC_RESET_CG (1<<7)
#define EIGS_ORTHOGONAL_PRECONDITIONER (1<<8)

/* default flags: what we think works best most of the time: */
#define EIGS_DEFAULT_FLAGS (EIGS_RESET_CG | EIGS_REORTHOGONALIZE)

typedef struct evectconstraint_chain_struct {
     evectconstraint C;
     void *constraint_data;
     struct evectconstraint_chain_struct *next;
} evectconstraint_chain;

extern evectconstraint_chain *evect_add_constraint(evectconstraint_chain 
						   *constraints,
						   evectconstraint C,
						   void *constraint_data);
extern void evect_destroy_constraints(evectconstraint_chain *constraints);
extern void evectconstraint_chain_func(evectmatrix X, void *data);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* EIGENSOLVER_H */
