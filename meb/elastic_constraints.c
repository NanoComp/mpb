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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../config.h"
#include <check.h>

#include <mpiglue.h>
#include "elastic.h"

/**************************************************************************/

/* function to call z and y parity constraints, if necessary */
void elastic_parity_constraint(evectmatrix X, void *data)
{
     elastic_data *d = (elastic_data *) data;

     CHECK(d, "null elastic data pointer!");
     CHECK(X.c == 3, "fields don't have 3 components!");

     CHECK(d->parity == NO_PARITY, "parity not yet supported in MEB");

/* TODO
     if (d->parity & (EVEN_Z_PARITY | ODD_Z_PARITY))
	  elastic_zparity_constraint(X, data);
     if (d->parity & (EVEN_Y_PARITY | ODD_Y_PARITY))
	  elastic_yparity_constraint(X, data);
*/
}

/**************************************************************************/

