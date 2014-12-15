noinst_LTLIBRARIES = libutil.la

libutil_la_SOURCES = check.h debug_malloc.c mpi_utils.c mpi_utils.h mpiglue.h sphere-quad.h

BUILT_SOURCES = sphere-quad.h

noinst_PROGRAMS = sphere_quad
sphere_quad_SOURCES = sphere-quad.c

if MAINTAINER_MODE

PRELUDE = "/* This file was automatically generated --- DO NOT EDIT */"
SPHERE_QUAD = $(top_builddir)/src/util/sphere_quad
sphere-quad.h: $(SPHERE_QUAD)$(EXEEXT)
	(echo $(PRELUDE); echo; $(SPHERE_QUAD)) > $@

endif
