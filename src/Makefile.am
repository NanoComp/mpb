SUBDIRS = util matrices matrixio maxwell

lib_LTLIBRARIES = libmpb@MPB_SUFFIX@.la

libmpb@MPB_SUFFIX@_la_SOURCES =
libmpb@MPB_SUFFIX@_la_LIBADD = util/libutil.la matrices/libmatrices.la maxwell/libmaxwell.la
libmpb@MPB_SUFFIX@_la_LDFLAGS = -no-undefined -version-info @SHARED_VERSION_INFO@

BUILT_SOURCES = mpb@MPB_SUFFIX@.h
include_HEADERS = mpb@MPB_SUFFIX@.h
pkginclude_HEADERS = matrices/eigensolver.h matrices/matrices.h	\
matrices/scalar.h maxwell/maxwell.h

mpb@MPB_SUFFIX@.h: mpbconf.h
	cp -f mpbconf.h $@

EXTRA_DIST = mpbconf.h.in
