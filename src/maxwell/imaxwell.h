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

/* internal header for maxwell.h, depends on config.h */

#ifndef IMAXWELL_H
#define IMAXWELL_H

#include "config.h"

#include "maxwell.h"

#if defined(HAVE_LIBFFTW)
#  include <fftw.h>
#  include <rfftw.h>
#  ifdef HAVE_MPI
#    include <fftw_mpi.h>
#    include <rfftw_mpi.h>
#  endif
#elif defined(HAVE_LIBDFFTW)
#  include <dfftw.h>
#  include <drfftw.h>
#  ifdef HAVE_MPI
#    include <dfftw_mpi.h>
#    include <drfftw_mpi.h>
#  endif
#elif defined(HAVE_LIBSFFTW)
#  include <sfftw.h>
#  include <srfftw.h>
#  ifdef HAVE_MPI
#    include <sfftw_mpi.h>
#    include <srfftw_mpi.h>
#  endif
#elif defined(HAVE_LIBXFFTW)
#  include <xfftw.h>
#  include <xrfftw.h>
#  ifdef HAVE_MPI
#    include <xfftw_mpi.h>
#    include <xrfftw_mpi.h>
#  endif
#endif

#if defined(HAVE_LIBFFTW) || defined(HAVE_LIBDFFTW) || defined(HAVE_LIBSFFTW) || defined(HAVE_LIBXFFTW)
#  define HAVE_FFTW 1
#endif

#ifdef HAVE_FFTW
#  ifdef HAVE_MPI
#    ifdef SCALAR_COMPLEX
     typedef fftwnd_mpi_plan fftplan;
#    else
     typedef rfftwnd_mpi_plan fftplan;
#    endif
#  else
#    ifdef SCALAR_COMPLEX
     typedef fftwnd_plan fftplan;
#    else
     typedef rfftwnd_plan fftplan;
#    endif
#  endif
#endif

#endif /* IMAXWELL_H */
