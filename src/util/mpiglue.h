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

#ifndef MPIGLUE_H
#define MPIGLUE_H

/* This header file stands between our code and MPI.  If MPI is present,
   we just #include <mpi.h>.  Otherwise, we use no-op stubs for 
   MPI routines. */

#ifdef HAVE_MPI

#include <mpi.h>

#else /* don't have MPI */

#include <check.h>

#define MPI_Allreduce(sb, rb, n, t, op, comm) \
CHECK((sb) == (rb), "MPI_Allreduce stub doesn't work for sendbuf != recvbuf")

#define MPI_Abort(comm, errcode) exit(errcode)

#define MPI_Barrier(comm) 0

#define MPI_Comm_rank(comm, rankp) *(rankp) = 0

#endif /* HAVE_MPI */

#endif /* MPIGLUE_H */
