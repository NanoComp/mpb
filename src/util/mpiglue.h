/* Copyright (C) 1999-2012, Massachusetts Institute of Technology.
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

#ifndef MPIGLUE_H
#define MPIGLUE_H

/* This header file stands between our code and MPI.  If MPI is present,
   we just #include <mpi.h>.  Otherwise, we use no-op stubs for 
   MPI routines. */

#ifdef HAVE_MPI

#include <mpi.h>

typedef double mpiglue_clock_t;
#define MPIGLUE_CLOCK MPI_Wtime()
#define MPIGLUE_CLOCK_DIFF(t2, t1) ((t2) - (t1))

#define mpi_allreduce(sb, rb, n, ctype, t, op, comm) \
     MPI_Allreduce(sb,rb,n,t,op,comm)

#else /* don't have MPI */

#include <string.h>
#include <time.h>
#include <check.h>

#define MPI_Init(argc,argv) 0
#define MPI_Finalize() 0

/* Stub for the MPI_Allreduce function, differing only by the addition
   of a "ctype" parameter that is the C type for the buffers.  This
   routine just loops over the buffers and assigns sb to rb, since
   all of the MPI reduce operations become simple assignments on a
   single processor.   (Our original implementation required rb == sb,
   but it seems that MPI implementations do not allow this.  The MPI 2.0
   specification supports an MPI_IN_PLACE constant that you can pass
   for sb in order to be in-place, but I don't want to require that. */
#define mpi_allreduce(sb, rb, n, ctype, t, op, comm) { \
     CHECK((sb) != (rb), "MPI_Allreduce doesn't work for sendbuf == recvbuf");\
     memcpy((rb), (sb), (n) * sizeof(ctype)); \
}

#define MPI_Bcast(b, n, t, root, comm) 0

#define MPI_Abort(comm, errcode) exit(errcode)

#define MPI_Barrier(comm) 0

#define MPI_Comm_rank(comm, rankp) *(rankp) = 0
#define MPI_Comm_size(comm, sizep) *(sizep) = 1

#define MPI_Send(sb,n,t, r,tag, comm) \
CHECK(0, "MPI_Send stub is non-functional");

#define MPI_Recv(sb,n,t, r,tag, comm,statp) \
CHECK(0, "MPI_Recv stub is non-functional");

typedef int mpiglue_status_t;
#define MPI_Status mpiglue_status_t

typedef clock_t mpiglue_clock_t;
#define MPIGLUE_CLOCK clock()
#define MPIGLUE_CLOCK_DIFF(t2, t1) (((t2) - (t1)) * 1.0 / CLOCKS_PER_SEC)

#endif /* HAVE_MPI */

#endif /* MPIGLUE_H */
