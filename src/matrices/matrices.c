/* Copyright (C) 1999 Massachusetts Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

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
