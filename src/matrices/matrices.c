#include <stdlib.h>
#include <stdio.h>

#include <config.h>
#include <check.h>

#include "matrices.h"

/* Basic operations: allocation, deallocation, etcetera. */

evectmatrix create_evectmatrix(int N, int c, int p,
			       int localN, int Nstart, int allocN)
{
     evectmatrix X;
 
     CHECK(localN <= N && allocN >= localN && Nstart <= localN,
	   "invalid N arguments");
    
     X.N = N;
     X.localN = localN;
     X.Nstart = Nstart;
     X.allocN = allocN;
     X.c = c;
     
     X.n = localN * c;
     X.p = p;
     
     if (allocN > 0) {
	  X.data = (scalar*) malloc(sizeof(scalar) * allocN * c * p);
	  CHECK(X.data, "out of memory");
     }
     else
	  X.data = NULL;

     return X;
}

void destroy_evectmatrix(evectmatrix X)
{
     free(X.data);
}

sqmatrix create_sqmatrix(int p)
{
     sqmatrix X;

     X.p = p;

     X.data = (scalar*) malloc(sizeof(scalar) * p * p);
     CHECK(X.data, "out of memory");

     return X;
}

void destroy_sqmatrix(sqmatrix X)
{
     free(X.data);
}
