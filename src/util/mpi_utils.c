/* Copyright (C) 1999, 2000, 2001 Massachusetts Institute of Technology.
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
#include <stdarg.h>

#include "../config.h"
#include <mpiglue.h>

#include "mpi_utils.h"

/* die when fatal errors occur */
void mpi_die(const char *template, ...)
{
     va_list ap;

     va_start(ap, template);
     vfprintf(stderr, template, ap);
     va_end(ap);

     MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}

void mpi_one_fprintf(FILE *f, const char *template, ...)
{
     if (mpi_is_master()) {
	  va_list ap;
	  va_start(ap, template);
	  vfprintf(f, template, ap);
	  va_end(ap);
     }
}

int mpi_is_master(void)
{
     int process_rank;
     MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
     return (process_rank == 0);
}
