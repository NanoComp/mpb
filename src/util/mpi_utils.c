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

#if !defined(_GNU_SOURCE)
#define _GNU_SOURCE
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include "config.h"
#include <check.h>
#include <mpiglue.h>

#include "mpi_utils.h"

#ifdef HAVE_MPI
MPI_Comm mpb_comm = MPI_COMM_WORLD;
#else
int mpb_comm = 0;
#endif

/* Simple, somewhat hackish API to allow user to run multiple simulations
   in parallel in the same MPI job.  The user calls

   mygroup = divide_parallel_processes(numgroups);

   to divide all of the MPI processes into numgroups equal groups,
   and to return the index (from 0 to numgroups-1) of the current group.
   From this point on, all fields etc. that you create and all
   calls using mpb_comm will only communicate within your group of
   processes.

   However, there are two calls that you can use to switch back to
   globally communication among all processes:

   begin_global_communications();
   ....do stuff....
   end_global_communications();

   It is important not to mix the two types; e.g. you cannot solve
   a field created in the local group in global mode, or vice versa.
*/

#ifdef HAVE_MPI
static MPI_Comm mpb_comm_save = MPI_COMM_WORLD;
#endif

void end_divide_parallel(void) {
#ifdef HAVE_MPI
  if (mpb_comm != MPI_COMM_WORLD) MPI_Comm_free(&mpb_comm);
  if (mpb_comm_save != MPI_COMM_WORLD) MPI_Comm_free(&mpb_comm_save);
  mpb_comm = mpb_comm_save = MPI_COMM_WORLD;
#endif
}

int divide_parallel_processes(int numgroups) {
#ifdef HAVE_MPI
  int sz, rank, mygroup;
  end_divide_parallel();
  MPI_Comm_size(mpb_comm, &sz);
  CHECK(numgroups > 0, "numgroups must be > 0");
  CHECK(numgroups <= sz, "tried to split into more groups than processes");
  MPI_Comm_rank(mpb_comm, &rank);
  mygroup = (rank * numgroups) / sz;
  MPI_Comm_split(MPI_COMM_WORLD, mygroup, rank, &mpb_comm);
  return mygroup;
#else
  CHECK(numgroups != 1, "tried to split into more groups than processes");
  return 0;
#endif
}

void begin_global_communications(void) {
#ifdef HAVE_MPI
  mpb_comm_save = mpb_comm;
  mpb_comm = MPI_COMM_WORLD;
#endif
}

void end_global_communications(void) {
#ifdef HAVE_MPI
  mpb_comm = mpb_comm_save;
  mpb_comm_save = MPI_COMM_WORLD;
#endif
}

int my_global_rank() {
#ifdef HAVE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
#else
  return 0;
#endif
}

/* die when fatal errors occur */
void mpi_die(const char *template, ...) {
  va_list ap;

  va_start(ap, template);
  vfprintf(stderr, template, ap);
  va_end(ap);

  MPI_Abort(mpb_comm, EXIT_FAILURE);
}

void (*mpb_printf_callback)(const char *s) = NULL;

/* Like printf, except only does anything on master process. */
void mpi_one_printf(const char *template, ...) {
  if (mpi_is_master()) {
    va_list ap;
    va_start(ap, template);
    if (mpb_printf_callback) {
      char *s;
      vasprintf(&s, template, ap);
      mpb_printf_callback(s);
      free(s);
    }
    else {
      vprintf(template, ap);
    }
    va_end(ap);
    fflush(stdout);
  }
}

/* Like fprintf, except only does anything on master process. */
void mpi_one_fprintf(FILE *f, const char *template, ...) {
  if (mpi_is_master()) {
    va_list ap;
    va_start(ap, template);
    vfprintf(f, template, ap);
    va_end(ap);
    fflush(f);
  }
}

/* Return whether we are the master process (rank == 0). */
int mpi_is_master(void) {
  int process_rank;
  MPI_Comm_rank(mpb_comm, &process_rank);
  return (process_rank == 0);
}

/* When debugging, checks to see that x is the same over all processes,
   and abort the program if it is not. */
void mpi_assert_equal(double x) {
#ifdef DEBUG
  double xmin, xmax;

  mpi_allreduce(&x, &xmin, 1, double, MPI_DOUBLE, MPI_MIN, mpb_comm);
  mpi_allreduce(&x, &xmax, 1, double, MPI_DOUBLE, MPI_MAX, mpb_comm);
  CHECK(xmin == x && xmax == x, "mpi_assert_equal failure");
#else
  (void)x; /* unused */
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

void mpi_begin_critical_section(int tag) {
  int process_rank;
  MPI_Comm_rank(mpb_comm, &process_rank);
  if (process_rank > 0) { /* wait for a message before continuing */
    MPI_Status status;
    int recv_tag = tag - 1; /* initialize to wrong value */
    MPI_Recv(&recv_tag, 1, MPI_INT, process_rank - 1, tag, mpb_comm, &status);
    CHECK(recv_tag == tag, "invalid tag received");
  }
}

void mpi_end_critical_section(int tag) {
  int process_rank, num_procs;
  MPI_Comm_rank(mpb_comm, &process_rank);
  MPI_Comm_size(mpb_comm, &num_procs);
  if (process_rank != num_procs - 1) { /* send a message to next process */
    MPI_Send(&tag, 1, MPI_INT, process_rank + 1, tag, mpb_comm);
  }
}
