if !MPI
bin_PROGRAMS = mpb@MPB_SUFFIX@-data
endif

mpb@MPB_SUFFIX@_data_SOURCES = mpb-data.c
mpb@MPB_SUFFIX@_data_LDADD = $(top_builddir)/src/matrixio/libmatrixio.a $(top_builddir)/src/libmpb@MPB_SUFFIX@.la
mpb@MPB_SUFFIX@_data_CPPFLAGS = -I$(top_srcdir)/src/util -I$(top_srcdir)/src/matrices -I$(top_srcdir)/src/matrixio

dist_man_MANS = mpb-data.1
