#ifndef MATRICES_H
#define MATRICES_H

#include "scalar.h"

typedef struct {
     int N, localN, Nstart, allocN;
     int c;
     int n, p;
     scalar *data;
} evectmatrix;

typedef struct {
     int p;
     scalar *data;
} sqmatrix;

/* general creation/destruction operations: */

extern evectmatrix create_evectmatrix(int N, int c, int p,
				      int localN, int Nstart, int allocN);
extern void destroy_evectmatrix(evectmatrix X);
extern sqmatrix create_sqmatrix(int p);
extern void destroy_sqmatrix(sqmatrix X);

/* evectmatrix operations, defined in evectmatrix.c: */

extern void evectmatrix_copy(evectmatrix X, evectmatrix Y);
extern void evectmatrix_aXpbY(real a, evectmatrix X, real b, evectmatrix Y);
extern void evectmatrix_XeYS(evectmatrix X, evectmatrix Y,
			     sqmatrix S, short sdagger);
extern void evectmatrix_XpaYS(evectmatrix X, real a, evectmatrix Y,sqmatrix S);
extern void evectmatrix_XtX(sqmatrix U, evectmatrix X);
extern void evectmatrix_XtY(sqmatrix U, evectmatrix X, evectmatrix Y);
extern scalar evectmatrix_traceXtY(evectmatrix X, evectmatrix Y);

/* sqmatrix operations, defined in sqmatrix.c: */

extern scalar sqmatrix_trace(sqmatrix U);
extern scalar sqmatrix_traceAtB(sqmatrix A, sqmatrix B);
extern void sqmatrix_AeBC(sqmatrix A, sqmatrix B, short bdagger,
			  sqmatrix C, short cdagger);
extern void sqmatrix_ApaB(sqmatrix A, real a, sqmatrix B);
extern void sqmatrix_invert(sqmatrix U, short compute_Uinv);
extern void sqmatrix_eigensolve(sqmatrix U, real *eigenvals, sqmatrix W);
extern void sqmatrix_sqrt(sqmatrix Usqrt, sqmatrix U, sqmatrix W);

#endif /* MATRICES_H */
