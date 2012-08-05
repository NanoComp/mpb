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

#ifndef MPI_UTILS_H
#define MPI_UTILS_H

#include <stdio.h>

extern void mpi_die(const char *template, ...)
#ifdef __GNUC__
     __attribute__ ((format (printf, 1, 2)))
#endif
;

extern void mpi_one_fprintf(FILE *f, const char *template, ...)
#ifdef __GNUC__
     __attribute__ ((format (printf, 2, 3)))
#endif
;

extern void mpi_one_printf(const char *template, ...)
#ifdef __GNUC__
     __attribute__ ((format (printf, 1, 2)))
#endif
;

extern int mpi_is_master(void);

extern void mpi_assert_equal(double x);

extern void mpi_begin_critical_section(int tag);
extern void mpi_end_critical_section(int tag);

/* "in-place" Allreduce wrapper for reducing a single value */
#define mpi_allreduce_1(b, ctype, t, op, comm) { \
     ctype bbbb = *(b); \
     mpi_allreduce(&bbbb, (b), 1, ctype, t, op, comm); \
}

#endif /* MPI_UTILS_H */
