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
#include <check.h>
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

/* Like printf, except only does anything on master process. */
void mpi_one_printf(const char *template, ...)
{
     if (mpi_is_master()) {
	  va_list ap;
	  va_start(ap, template);
	  vprintf(template, ap);
	  va_end(ap);
     }
}

/* Like fprintf, except only does anything on master process. */
void mpi_one_fprintf(FILE *f, const char *template, ...)
{
     if (mpi_is_master()) {
	  va_list ap;
	  va_start(ap, template);
	  vfprintf(f, template, ap);
	  va_end(ap);
     }
}

/* Return whether we are the master process (rank == 0). */
int mpi_is_master(void)
{
     int process_rank;
     MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
     return (process_rank == 0);
}

/* When debugging, checks to see that x is the same over all processes,
   and abort the program if it is not. */
void mpi_assert_equal(double x)
{
#ifdef DEBUG
     double xmin, xmax;

     mpi_allreduce(&x, &xmin, 1, double, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
     mpi_allreduce(&x, &xmax, 1, double, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
     CHECK(xmin == x && xmax == x, "mpi_assert_equal failure");
#endif
}

/* The following functions bracket a "critical section," a region
   of code that should be executed by only one process at a time.

   They work by having each process wait for a message from the
   previous process before starting. 

   Each critical section is passed an integer "tag"...ideally, this
   should be a unique identifier for each critical section so that
   messages from different critical sections don't get mixed up
   somehow. */

void mpi_begin_critical_section(int tag)
{
     int process_rank;
     MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
     if (process_rank > 0) { /* wait for a message before continuing */
	  MPI_Status status;
	  int recv_tag = tag - 1; /* initialize to wrong value */
	  MPI_Recv(&recv_tag, 1, MPI_INT, process_rank - 1, tag, 
		   MPI_COMM_WORLD, &status);
	  CHECK(recv_tag == tag, "invalid tag received");
     }
}

void mpi_end_critical_section(int tag)
{
     int process_rank, num_procs;
     MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
     MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
     if (process_rank != num_procs - 1) { /* send a message to next process */
	  MPI_Send(&tag, 1, MPI_INT, process_rank + 1, tag, 
		   MPI_COMM_WORLD);
     }
}
