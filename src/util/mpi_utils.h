#ifndef MPI_UTIL_H
#define MPI_UTIL_H

#include <stdio.h>

extern void mpi_die(const char *template, ...)
#ifdef __GNUC__
     __attribute__ ((format (printf, 1, 2)));
#endif

extern void mpi_one_fprintf(FILE *f, const char *template, ...)
#ifdef __GNUC__
     __attribute__ ((format (printf, 2, 3)));
#endif

#endif /* MPI_UTIL_H */
