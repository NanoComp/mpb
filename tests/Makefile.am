noinst_PROGRAMS = malloctest blastest eigs_test maxwell_test normal_vectors
EXTRA_DIST = blastest.real.out blastest.complex.out

LIBMPB = $(top_builddir)/src/libmpb@MPB_SUFFIX@.la
AM_CPPFLAGS = -I$(top_srcdir)/src/util -I$(top_srcdir)/src/matrices -I$(top_srcdir)/src/maxwell

malloctest_SOURCES = malloctest.c
malloctest_LDADD = $(LIBMPB)

blastest_SOURCES = blastest.c
blastest_LDADD = $(LIBMPB)

eigs_test_SOURCES = eigs_test.c
eigs_test_LDADD = $(LIBMPB)

maxwell_test_SOURCES = maxwell_test.c
maxwell_test_LDADD = $(LIBMPB)

normal_vectors_SOURCES = normal_vectors.c
normal_vectors_LDADD = $(LIBMPB)

check-local: blastest maxwell_test # malloctest eigs_test normal_vectors
	rm -f blastest.out.test
	./blastest | sed 's/\-0\.000\([ ,)]\)/ 0.000\1/g' | sed 's/\-0\.000$$/ 0.000/g' > blastest.out.test
	diff $(srcdir)/blastest.@SCALAR_TYPE@.out blastest.out.test
	./maxwell_test -1 -c 1e-9 -x 256 -E 1e-3
	@echo "**********************************************************"
	@echo "                       PASSED tests."
	@echo "**********************************************************"
