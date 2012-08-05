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

#ifndef CHECK_H
#define CHECK_H

#ifdef CHECK_DISABLE

#define CHECK(condition, message) /* do nothing */

#else /* not CHECK_DISABLE */

#include <mpi_utils.h>

#ifdef DEBUG
extern void check_breakpoint(void);
#define CHECK_BREAKPOINT check_breakpoint() /* defined in debug_malloc.c */
#else
#define CHECK_BREAKPOINT /* nothing */
#endif

#define CHECK(condition, message) do { \
     if (!(condition))  { \
          CHECK_BREAKPOINT; \
          mpi_die("CHECK failure on line %d of " __FILE__ ": " \
		  message "\n", __LINE__); \
     } \
} while (0)

#endif /* not CHECK_DISABLE */

#ifdef DEBUG_MALLOC
extern void *debug_malloc(size_t n);
extern void debug_free(void *p);
#  define malloc debug_malloc
#  define free debug_free
#endif

/* Macro to check whether a floating-point number contains a ridiculous
   value.  Note that x != x if and only if x is a NaN. */
#define BADNUM(x) ((x) != (x) || (x) > 1e50 || (x) < -1e50)

extern void debug_output_malloc_count(void);
extern void debug_check_memory_leaks(void);

#define CHK_MALLOC(p, t, n) {                                         \
     size_t CHK_MALLOC_n_tmp = (n);                                   \
     (p) = (t *) malloc(sizeof(t) * CHK_MALLOC_n_tmp);                \
     CHECK((p) || CHK_MALLOC_n_tmp == 0, "out of memory!");           \
}

#endif /* CHECK_H */
