noinst_LTLIBRARIES = libmatrices.la
EXTRA_DIST = README

libmatrices_la_SOURCES = blasglue.c blasglue.h eigensolver.c		\
eigensolver.h eigensolver_davidson.c eigensolver_utils.c evectmatrix.c	\
linmin.c linmin.h matrices.c matrices.h minpack2-linmin.c scalar.h	\
sqmatrix.c
libmatrices_la_CPPFLAGS = -I$(srcdir)/../util
