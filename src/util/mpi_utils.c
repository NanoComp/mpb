#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include <config.h>
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
     va_list ap;
     int process_rank;

     MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
     
     if (process_rank == 0) {
	  va_start(ap, template);
	  vfprintf(f, template, ap);
	  va_end(ap);
     }
}
