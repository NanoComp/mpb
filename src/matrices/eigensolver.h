#ifndef EIGENSOLVER_H
#define EIGENSOLVER_H

#include <matrices.h>

typedef void (*evectoperator) (evectmatrix Xin, evectmatrix Xout,
			       void *data,
			       int is_current_eigenvector);

typedef void (*evectpreconditioner) (evectmatrix Xin, evectmatrix Xout,
				     void *data,
				     evectmatrix Y, real *eigenvals);

typedef void (*evectpreconditioner_data_updater) (void *data,
						  sqmatrix Yrotation);

extern void eigensolver(evectmatrix Y, real *eigenvals,
			evectoperator A, void *Adata,
			evectpreconditioner C, void *Cdata,
			evectpreconditioner_data_updater Cdata_update,
			evectmatrix Work[], int nWork,
			real tolerance, int *num_iterations);

#endif /* EIGENSOLVER_H */
