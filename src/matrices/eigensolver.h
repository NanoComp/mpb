#ifndef EIGENSOLVER_H
#define EIGENSOLVER_H

#include <matrices.h>

typedef void (*evectoperator) (evectmatrix Xin, evectmatrix Xout,
			       void *data);

extern void eigensolver(evectmatrix Y, real *eigenvals,
			evectoperator A, void *Adata,
			evectoperator C, void *Cdata,
			evectmatrix Work[], int nWork,
			real tolerance, int *num_iterations);

#endif /* EIGENSOLVER_H */
