noinst_PROGRAMS = malloctest blastest eigs_test maxwell_test
if WITH_LIBCTLGEOM
noinst_PROGRAMS += normal_vectors
endif
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
normal_vectors_LDADD = -lctlgeom $(LIBMPB)
normal_vectors_CPPFLAGS = $(CTLGEOM_H_CPPFLAG) $(AM_CPPFLAGS)

blastest.out: blastest
	./blastest | sed 's/\-0\.000\([ ,)]\)/ 0.000\1/g' | sed 's/\-0\.000$$/ 0.000/g' > $@
	diff $(srcdir)/blastest.@SCALAR_TYPE@.out $@

maxwell_test.out: maxwell_test
	./maxwell_test -1 -c 1e-9 -x 256 -E 1e-3 > $@

maxwell_test_2.out: maxwell_test
	./maxwell_test -1 -c 1e-9 -x 256 -E 1e-3 -e -k 0.4 -n 1 > $@

if !MPI
MAXWELL_TEST_OUT=maxwell_test.out maxwell_test_2.out
endif

check-local: blastest.out $(MAXWELL_TEST_OUT)
	@echo "**********************************************************"
	@echo "                       PASSED tests."
	@echo "**********************************************************"

clean-local:
	rm -f blastest.out $(MAXWELL_TEST_OUT)
