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
