SUBDIRS = util matrices matrixio maxwell

lib_LTLIBRARIES = libmpb@MPB_SUFFIX@.la

libmpb@MPB_SUFFIX@_la_SOURCES =
libmpb@MPB_SUFFIX@_la_LIBADD = util/libutil.la matrices/libmatrices.la maxwell/libmaxwell.la
libmpb@MPB_SUFFIX@_la_LDFLAGS = -no-undefined -version-info @SHARED_VERSION_INFO@

pkginclude_HEADERS = matrices/eigensolver.h matrices/matrices.h	\
matrices/scalar.h maxwell/maxwell.h
