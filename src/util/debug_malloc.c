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

/* Debugging versions of malloc & free, to help catch common errors.
   The routines in this file were adapted from similar routines in the
   FFTW library (www.fftw.org), written by Matteo Frigo and Steven
   G. Johnson.  FFTW is also GPLed and is also copyrighted by MIT.  */

#include <stdio.h>
#include <stdlib.h>

#include "../config.h"

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#include "mpi_utils.h"
#include "check.h"

/**********************************************************
 *   DEBUGGING CODE
 **********************************************************/

#ifdef DEBUG

void check_breakpoint(void)
{
     /* this function is only here so that we can drop into a breakpoint
	in a debugger when a CHECK macro fails... */
}

#endif

#ifdef DEBUG_MALLOC

#  undef malloc
#  undef free

/*
 * debugging malloc/free.  Initialize every malloced and freed area to
 * random values, just to make sure we are not using uninitialized
 * pointers.  Also check for writes past the ends of allocated blocks,
 * and a couple of other things.
 *
 */

static int debug_malloc_cnt = 0;
static int debug_malloc_total = 0;

#define MAGIC 0xABadCafe
#define MMAGIC (((int) MAGIC) < 0 ? ((int) MAGIC) : -((int) MAGIC))
#define PAD_FACTOR 2
#define TWOINTS (2 * sizeof(int))

#define VERBOSE_ALLOCATION 0

#if VERBOSE_ALLOCATION
#define WHEN_VERBOSE(a) a
#else
#define WHEN_VERBOSE(a)
#endif

void *debug_malloc(size_t n)
{
     char *p;
     int i;

     WHEN_VERBOSE(mpi_one_fprintf(stdout,"DEBUG_MALLOC %d\n", n));

     if (n == 0)
	  mpi_one_fprintf(stderr, "(Allocating a block of zero size.)\n");

     debug_malloc_total += n;

     p = (char *) malloc(PAD_FACTOR * n + TWOINTS);
     CHECK(p, "debug_malloc: out of memory\n");

     /* store the size in a known position */
     ((int *) p)[0] = n;
     ((int *) p)[1] = MAGIC;
     for (i = 0; i < PAD_FACTOR * n; ++i)
	  p[i + TWOINTS] = (char) (i ^ 0xDEADBEEF);

     ++debug_malloc_cnt;

     /* skip the size we stored previously */
     return (void *) (p + TWOINTS);
}

void debug_free(void *p)
{
     char *q = ((char *) p) - TWOINTS;

     CHECK(p, "debug_free: tried to free NULL pointer!\n");
     CHECK(q, "debug_free: tried to free NULL+TWOINTS pointer!\n");

     {
	  int n = ((int *) q)[0];
	  int magic = ((int *) q)[1];
	  int i;

	  WHEN_VERBOSE(mpi_one_fprintf(stdout,"DEBUG_FREE %d\n", n));
	  CHECK(n != MMAGIC, "Tried to free a freed pointer!\n");
	  *((int *) q) = MMAGIC; /* to detect duplicate free's */

	  CHECK(magic == MAGIC, "Wrong magic in debug_free()!\n");
	  ((int *) q)[1] = ~MAGIC;

	  CHECK(n >= 0, "Tried to free block with corrupt size descriptor!\n");

	  debug_malloc_total -= n;

	  CHECK(debug_malloc_total >= 0, 
		"debug_malloc_total went negative!\n");

	  /* check for writing past end of array: */
	  for (i = n; i < PAD_FACTOR * n; ++i)
	       if (q[i + TWOINTS] != (char) (i ^ 0xDEADBEEF))
		    mpi_die("Byte %d past end of array has changed!\n"
			    "Array bounds overwritten!\n",
			    i - n + 1);
	  for (i = 0; i < PAD_FACTOR * n; ++i)
	       q[i + TWOINTS] = (char) (i ^ 0xBEEFDEAD);

	  --debug_malloc_cnt;
	  free(q);
     }
}

#endif /* DEBUG */

/* output current memory usage: */
void debug_output_malloc_count(void)
{
#ifdef DEBUG_MALLOC
     mpi_one_fprintf(stderr, "malloc: %d blocks, %g kB\n", 
		     debug_malloc_cnt, debug_malloc_total / 1024.0);
#endif
}

/* check for memory leaks when debugging */
void debug_check_memory_leaks(void)
{
#ifdef DEBUG_MALLOC
     if (debug_malloc_cnt || debug_malloc_total)
	  mpi_die("MEMORY LEAK!!!\n"
		  "number of unbalanced malloc calls = %d\n"
		  "total leaked bytes = %d\n",
		  debug_malloc_cnt,
		  debug_malloc_total);
#endif
}
