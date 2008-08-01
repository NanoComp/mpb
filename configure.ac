# Process this file with autoconf to produce a configure script.
AC_INIT(mpb, 1.5, stevenj@alum.mit.edu)
AC_CONFIG_SRCDIR([src/matrices/matrices.c])
AC_CONFIG_HEADER([config.h src/mpbconf.h])

# Shared-library version number; indicates api compatibility, and is
# not the same as the "public" version number.  (Don't worry about this
# except for public releases.)
SHARED_VERSION_INFO="0:0:0"

AM_INIT_AUTOMAKE(1.6)
AC_SUBST(SHARED_VERSION_INFO)
AM_ENABLE_SHARED(no) dnl shared libs cause too many headaches to be default
AC_PROG_LIBTOOL

###########################################################################
MPB_VERSION=$PACKAGE_VERSION
AC_SUBST(MPB_VERSION)

AC_MSG_CHECKING([for vendor f77 to be used instead of g77])
AC_MSG_RESULT()
AC_CHECK_PROGS(F77, f77 xlf77 cf77 fort77 f90 xlf90)
AC_PROG_F77

# Get C compiler.
AC_MSG_CHECKING([for vendor cc to be used instead of gcc])
AC_MSG_RESULT()
AC_CHECK_PROGS(CC, cc xlc) # prefer vendor cc, to stay in sync with Fortran
AC_PROG_CC

AC_ARG_ENABLE(single, [AC_HELP_STRING([--enable-single],[compile for single precision])], enable_single=$enableval, enable_single=no)
if test "$enable_single" = "yes"; then
        AC_DEFINE(SCALAR_SINGLE_PREC,1,[Define to use single precision])
	echo "*********************** SINGLE PRECISION ***********************"
fi

AC_ARG_ENABLE(long-double, [AC_HELP_STRING([--enable-long-double],[compile for long-double precision])], enable_long_double=$enableval, enable_long_double=no)
if test "$enable_long_double" = "yes"; then
        AC_DEFINE(SCALAR_LONG_DOUBLE_PREC,1,[Define to use long-double precision])
	echo "****************** LONG-DOUBLE PRECISION ***********************"
fi

AC_ARG_ENABLE(debug, [AC_HELP_STRING([--enable-debug],[compile for debugging])], ok=$enableval, ok=no)
if test "$ok" = "yes"; then
	CFLAGS="-g"
	AC_DEFINE(DEBUG,1,[Define to turn on debugging checks])
fi

AC_ARG_ENABLE(debug-malloc, [AC_HELP_STRING([--enable-debug-malloc],[use debug malloc/free])], ok=$enableval, ok=no)
if test "$ok" = "yes"; then
	AC_DEFINE(DEBUG_MALLOC,1,[Define to use debugging malloc/free])
fi

AC_ARG_ENABLE(prof, [AC_HELP_STRING([--enable-prof],[compile for profiling])], ok=$enableval, ok=no)
if test "$ok" = "yes"; then
	CFLAGS="$CFLAGS -pg"
	AC_DEFINE(ENABLE_PROF,1,[Define when using the profiler tool])
fi

AC_ARG_ENABLE(checks, [AC_HELP_STRING([--disable-checks],[disable runtime checks])], ok=$enableval, ok=yes)
if test "$ok" = "no"; then
	AC_DEFINE(CHECK_DISABLE,1,[Define to disable sanity checks in code])
fi

AC_ARG_ENABLE(kottke, [AC_HELP_STRING([--disable-kottke],[disable new anisotropic smoothing based on Kottke algorithm])], ok=$enableval, ok=yes)
if test "$ok" = "yes"; then
	AC_DEFINE(KOTTKE,1,[Define to enable Kottke anisotropic smoothing])
fi

AC_ARG_WITH(inv-symmetry,[AC_HELP_STRING([--with-inv-symmetry],[take advantage of (and require) inv. sym.])], inv_sym=$withval, inv_sym=no)
if test "$inv_sym" = "no"; then
        AC_DEFINE(SCALAR_COMPLEX,1,[Define to use complex fields and
	                            not to require inversion symmetry])
	SCALAR_TYPE=complex
else
	SCALAR_TYPE=real
	# Make installed program mpbi instead of mpb, but don't change mpb-data
	MPB_SUFFIX="i"
	echo "*********************** INVERSION SYM. ***********************"
fi
AC_SUBST(SCALAR_TYPE)

AC_ARG_WITH(hermitian-eps,[AC_HELP_STRING([--with-hermitian-eps],[allow complex-Hermitian dielectric tensors])], hermitian_epsilon=$withval, hermitian_epsilon=no)
if test "$hermitian_epsilon" = "yes"; then
        AC_DEFINE(WITH_HERMITIAN_EPSILON,1,[Define to support Hermitian/complex
                                            dielectric tensors.])
	echo "*********************** HERMITIAN EPS. ***********************"
fi

# Add lots of compiler warnings to check for if we are using gcc:
# (The variable $GCC is set to "yes" by AC_PROG_CC if we are using gcc.)
if test "$GCC" = "yes"; then
	CFLAGS="$CFLAGS -Wall -W -Wbad-function-cast -Wcast-qual -Wpointer-arith -Wcast-align -pedantic"
fi

##############################################################################
#                       Checks for libraries:
##############################################################################

AC_ARG_WITH(efence, [AC_HELP_STRING([--with-efence],[use ElectricFence library])], ok=$withval, ok=no)
if test "$ok" = "yes"; then
	AC_CHECK_LIB(efence, malloc)
fi

##############################################################################
AC_CHECK_LIB(m, sqrt)

# If we have a Fortran compiler, link to its libraries, since these
# are probably needed to link with the BLAS/LAPACK libraries.
if test x != x"$F77"; then
	AC_F77_LIBRARY_LDFLAGS
	AC_F77_WRAPPERS
fi

##############################################################################
# In principle, we should abort if some of the following libraries are
# missing, or perhaps only compile some of the subdirectories.  For now,
# though, we do nothing (and wait for the compilation errors).

##############################################################################
AC_ARG_WITH(mpi, [AC_HELP_STRING([--with-mpi],[enable MPI parallelization])], with_mpi=$withval, with_mpi=no)
AC_ARG_WITH(fftw2, [AC_HELP_STRING([--with-fftw2],[prefer FFTW2 to FFTW3])], with_fftw2=$withval, with_fftw2=no)

# Check for FFTW libraries
fftw3_func=execute
test "x$with_mpi" = xno || fftw3_func=notarealfunction # skip FFTW3 with MPI
test "x$with_fftw2" = xyes && fftw3_func=notarealfunction # prefer FFTW2
if test "$enable_single" = "yes"; then
        AC_CHECK_LIB(fftw3f, fftwf_$fftw3_func, [], [
	AC_CHECK_LIB(sfftw, fftw)
	if test x"$inv_sym" = xyes; then
		AC_CHECK_LIB(srfftw, rfftw)
	fi])
elif test "$enable_long_double" = "yes"; then
        AC_CHECK_LIB(fftw3l, fftwl_$fftw3_func, [], [
	AC_CHECK_LIB(xfftw, fftw)
	if test x"$inv_sym" = xyes; then
		AC_CHECK_LIB(xrfftw, rfftw)
	fi])
else
        AC_CHECK_LIB(fftw3, fftw_$fftw3_func, [], [
	AC_CHECK_LIB(dfftw, fftw)
	if test x"$inv_sym" = xyes; then
		AC_CHECK_LIB(drfftw, rfftw)
	fi])
fi

if test x = x"`echo $LIBS | egrep 'l[[dsx]]fftw|fftw3'`"; then
	AC_CHECK_LIB(fftw, fftw)
	if test x"$inv_sym" = xyes; then
		AC_CHECK_LIB(rfftw, rfftw)
	fi
fi

if test x = x"`echo $LIBS | egrep 'l[[dsx]]*fftw'`"; then
	AC_MSG_ERROR([The FFTW libraries could not be found.  Make sure FFTW is installed, and that LDFLAGS and CPPFLAGS are set appropriately if necessary.])
fi

if test x"$inv_sym" = xyes; then
	if test x = x"`echo $LIBS | egrep 'l[[dsx]]*rfftw|lfftw3'`"; then
		AC_MSG_ERROR([The RFFTW libraries, which are required to compile MPB with inversion symmetry, could not be found.  These should have been installed as a part of FFTW.])
	fi
fi

##############################################################################
# Checks for BLAS and LAPACK libraries:

ACX_BLAS([], [AC_MSG_ERROR([BLAS was not found!])])
ACX_LAPACK([], [AC_MSG_ERROR([LAPACK was not found!])])

LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

##############################################################################
# Check for MPI library

if test "x$with_mpi" = "xyes"; then
	ACX_MPI([], [AC_MSG_ERROR(could not find mpi library for --with-mpi)])
        CC="$MPICC"
	LIBS="$MPILIBS $LIBS"

	# Make installed program mpb-mpi instead of mpb
	MPB_SUFFIX="${MPB_SUFFIX}_mpi"

	# Check for fftw_mpi library.
	if test "$enable_single" = "yes"; then
        	AC_CHECK_LIB(sfftw_mpi, fftwnd_mpi)
	        if test x"$inv_sym" = xyes; then
	                AC_CHECK_LIB(srfftw_mpi, rfftwnd_mpi)
	        fi
	elif test "$enable_long_double" = "yes"; then
        	AC_CHECK_LIB(xfftw_mpi, fftwnd_mpi)
	        if test x"$inv_sym" = xyes; then
	                AC_CHECK_LIB(xrfftw_mpi, rfftwnd_mpi)
	        fi
	else
	        AC_CHECK_LIB(dfftw_mpi, fftwnd_mpi)
	        if test x"$inv_sym" = xyes; then
	                AC_CHECK_LIB(drfftw_mpi, rfftwnd_mpi)
	        fi
	fi
	if test x = x"`echo $LIBS | egrep 'l[[ds]]fftw_mpi'`"; then
       		AC_CHECK_LIB(fftw_mpi, fftwnd_mpi)
	        if test x"$inv_sym" = xyes; then
	                AC_CHECK_LIB(rfftw_mpi, rfftwnd_mpi)
	        fi
	fi
	if test x = x"`echo $LIBS | egrep 'l[[ds]]*fftw_mpi'`"; then
		AC_MSG_ERROR([couldn't find fftw_mpi library for --with-mpi])
	fi
        if test x"$inv_sym" = xyes; then
		if test x = x"`echo $LIBS | egrep 'l[[ds]]*rfftw_mpi'`"; then
			AC_MSG_ERROR([couldn't find rfftw_mpi library for --with-mpi --with-inv-symmetry])
		fi
	fi

	AC_DEFINE(HAVE_MPI,1,[Define if you have & link an MPI library.])
fi

AM_CONDITIONAL(MPI, test "x$with_mpi" = "xyes")

##############################################################################
# Check for HDF5 library
AC_ARG_WITH(hdf5, [AC_HELP_STRING([--without-hdf5],[do not use HDF5 output])], ok=$withval, ok=yes)
if test "$ok" = "yes"; then
	AC_CHECK_LIB(z,deflate, [],[AC_MSG_WARN([zlib is required for HDF5!])])

	save_LIBS_0="$LIBS" # save, to check later if we found any library
	AC_CHECK_LIB(hdf5, H5Pcreate, [
		AC_CHECK_HEADERS(hdf5.h, [LIBS="-lhdf5 $LIBS"
                     AC_DEFINE(HAVE_HDF5,1,[Define if we have & link HDF5])])])
	if test x"$save_LIBS_0" = x"$LIBS"; then
		AC_MSG_WARN([Couldn't find the HDF5 library!!  Switching to --without-hdf5.])
	fi

        if test x != x"$MPILIBS"; then
		AC_CHECK_FUNCS(H5Pset_mpi H5Pset_fapl_mpio)
	fi
fi

##############################################################################
# Check for Guile library and flags:

AC_CHECK_PROG(guile_ok, guile, yes, no)
if test x"$guile_ok" = xno; then
	AC_MSG_ERROR([could not find guile program; check your PATH])
fi

AC_CHECK_PROG(HAVE_GUILE_CONFIG, guile-config, yes, no)
if test "$HAVE_GUILE_CONFIG" = "yes"; then
        CPPFLAGS="$CPPFLAGS `guile-config compile`"
        LIBS="`guile-config link` $LIBS"
	AC_MSG_CHECKING([if linking to guile works])
	AC_TRY_LINK_FUNC(gh_enter, AC_MSG_RESULT(yes), [
	AC_MSG_RESULT(no)
	AC_MSG_ERROR([guile-config is broken or Guile is not working.])
	])
else
        AC_CHECK_LIB(readline, readline)
        AC_CHECK_LIB(dl, dlopen)
        AC_CHECK_LIB(guile-ltdl, scm_lt_dlinit, [], [AC_CHECK_LIB(ltdl, dlopen)])
        AC_CHECK_LIB(guile, gh_eval_str, [], [AC_MSG_ERROR([linking to guile failed])])
fi

AC_CHECK_FUNCS(scm_make_smob_type scm_array_get_handle)
AC_MSG_CHECKING([for SCM_SMOB_PREDICATE])
AC_TRY_LINK([#include <guile/gh.h>
#include <libguile/smob.h>],
                    [int x; SCM_SMOB_PREDICATE(1,x);],
                    ok=yes, ok=no)
AC_MSG_RESULT($ok)
if test $ok = yes; then
	AC_DEFINE(HAVE_SCM_SMOB_PREDICATE, 1, [define if we have SCM_SMOB_PREDICATE])
fi
AC_MSG_CHECKING([for SCM_SMOB_DATA])
AC_TRY_LINK([#include <guile/gh.h>
#include <libguile/smob.h>],
                    [int x; SCM_SMOB_DATA(x);],
                    ok=yes, ok=no)
AC_MSG_RESULT($ok)
if test $ok = yes; then
	AC_DEFINE(HAVE_SCM_SMOB_DATA, 1, [define if we have SCM_SMOB_DATA])
fi
AC_MSG_CHECKING([for SCM_NEWSMOB])
AC_TRY_LINK([#include <guile/gh.h>
#include <libguile/smob.h>],
                    [int x; SCM_NEWSMOB(x,1,0);],
                    ok=yes, ok=no)
AC_MSG_RESULT($ok)
if test $ok = yes; then
	AC_DEFINE(HAVE_SCM_NEWSMOB, 1, [define if we have SCM_NEWSMOB])
fi

# Check how to activate Guile readline support:

ACTIVATE_READLINE=""
if test "x$cross_compiling" = xno; then
AC_MSG_CHECKING(how to activate readline in Guile)
ractivate="(use-modules (readline-activator)) (activate-readline)"
if guile -c "$ractivate" > /dev/null 2>&1; then
	AC_MSG_RESULT(readline-activator)
	ACTIVATE_READLINE="$ractivate"
else
	ractivate="(use-modules (ice-9 readline)) (activate-readline)"
	if guile -c "$ractivate" > /dev/null 2>&1; then
	        AC_MSG_RESULT(ice-9 readline)
	        ACTIVATE_READLINE="$ractivate"
	else
		AC_MSG_RESULT(cannot)
		ACTIVATE_READLINE=""
	fi
fi
fi
AC_SUBST(ACTIVATE_READLINE)

##############################################################################
# Check for libctl library and files

AC_ARG_WITH(libctl, [AC_HELP_STRING([--with-libctl=<dir>],[specify libctl directory])], ok=$withval, ok=yes)
if test "$ok" = "no"; then
        AC_MSG_ERROR([libctl is required and you may not disable it])
elif test "$ok" != "yes"; then
	LIBCTL_DIR="$ok"
fi

AC_MSG_CHECKING([for libctl dir])
if test x != x"$LIBCTL_DIR" -a ! -r "$LIBCTL_DIR/base/ctl.scm"; then
	LIBCTL_DIR=""
fi
if test x = x"$LIBCTL_DIR" -a -r /usr/local/share/libctl/base/ctl.scm; then
	LIBCTL_DIR="/usr/local/share/libctl"
fi
if test x = x"$LIBCTL_DIR" -a -r /usr/share/libctl3/base/ctl.scm; then
        LIBCTL_DIR="/usr/share/libctl3"
fi
if test x = x"$LIBCTL_DIR" -a -r /usr/share/libctl/base/ctl.scm; then
        LIBCTL_DIR="/usr/share/libctl"
fi
if test x = x"$LIBCTL_DIR"; then
	AC_MSG_ERROR([could not find libctl files; use --with-libctl=<dir>])
fi
AC_MSG_RESULT($LIBCTL_DIR)
AC_SUBST(LIBCTL_DIR)

AC_CHECK_PROG(gen_ctl_io_ok, gen-ctl-io, yes, no)
if test x"$gen_ctl_io_ok" = xno; then
	AC_MSG_ERROR([could not find gen-ctl-io program; check your PATH])
fi

GEN_CTL_IO=gen-ctl-io
AC_SUBST(GEN_CTL_IO)

# check for -lctl:
AC_CHECK_LIB(ctl, ctl_get_vector3, [], [AC_MSG_ERROR([Couldn't find the required libctl library.])])

# check for ctl.h (Fedora puts it in /usr/include/ctl, grrr)
CTL_H_CPPFLAG=""
AC_CHECK_HEADER(ctl.h, [CTL_H_CPPFLAG=" "], [
  save_CPPFLAGS=$CPPFLAGS
  AC_MSG_CHECKING([for alternate ctl.h location])
  for ctl_h_path in /usr/include/ctl /usr/local/include/ctl; do
     CPPFLAGS="$save_CPPFLAGS -I$ctl_h_path"
     AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[#include <ctl.h>]],[]), [CTL_H_CPPFLAG=-I$ctl_h_path; break])
  done
  AC_MSG_RESULT([$CTL_H_CPPFLAG])
  CPPFLAGS=$save_CPPFLAGS])
if test "x$CTL_H_CPPFLAG" = "x"; then
  AC_MSG_ERROR([Couldn't find the ctl.h header file for libctl.])
fi
AC_SUBST(CTL_H_CPPFLAG)

# Check libctl version >= LIBCTL_MAJOR.LIBCTL_MINOR.LIBCTL_BUGFIX
LIBCTL_MAJOR=3; LIBCTL_MINOR=0; LIBCTL_BUGFIX=0
AC_MSG_CHECKING([whether libctl version is at least ${LIBCTL_MAJOR}.${LIBCTL_MINOR}.${LIBCTL_BUGFIX}])
AC_EGREP_CPP(yes, [[
#include <ctl.h>
#if LIBCTL_MAJOR_VERSION > $LIBCTL_MAJOR || (LIBCTL_MAJOR_VERSION == $LIBCTL_MAJOR && (LIBCTL_MINOR_VERSION > $LIBCTL_MINOR || (LIBCTL_MINOR_VERSION == $LIBCTL_MINOR && LIBCTL_BUGFIX_VERSION >= $LIBCTL_BUGFIX)))
  yes
#endif
]], [AC_MSG_RESULT(ok)], [AC_MSG_ERROR([libctl version ${LIBCTL_MAJOR}.${LIBCTL_MINOR}.${LIBCTL_BUGFIX} or later is required])])

# On IRIX, basename/dirname functions in libctl/main.c require -lgen
AC_CHECK_LIB(gen, basename)

##############################################################################
# check for nlopt library and libctl support

with_nlopt=no
if test -r $LIBCTL_DIR/utils/nlopt.c -a -r $LIBCTL_DIR/utils/nlopt-constants.scm; then
   AC_CHECK_LIB(nlopt, nlopt_minimize, [with_nlopt=yes
					AC_DEFINE(HAVE_NLOPT,[1],[Have nlopt lib])
				        NLOPT_CPPFLAG=-DHAVE_NLOPT])
fi
AC_SUBST(NLOPT_CPPFLAG)
AM_CONDITIONAL(NLOPT, test "x$with_nlopt" = "xyes")

##############################################################################
# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS(unistd.h getopt.h nlopt.h)

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE

# Checks for library functions.
AC_CHECK_FUNCS(getopt strncmp)

##############################################################################
# Check to see if calling Fortran functions (in particular, the BLAS
# complex dot product functions) works from C:

AC_MSG_CHECKING([whether calling BLAS zdotc works])
AC_F77_FUNC(zdotc)
AC_TRY_RUN([
#include <stdio.h>
#include <stdlib.h>
typedef struct { double re, im; } cmplx;
#ifdef __cplusplus
extern "C"
#endif
extern cmplx $zdotc(int *, cmplx *, int *, cmplx *, int *);
int main(void) {
     cmplx x[2] = { {1,0}, {0,1} }, y[2] = { {1,0}, {1,0} };
     int n = 2, inc1 = 1, inc2 = 1;
     cmplx result = { 0, 0 };
     result = $zdotc(&n, x, &inc1, y, &inc2);
     if (result.re != 1 || result.im != -1) return EXIT_FAILURE;
     else return EXIT_SUCCESS;
}
], ok="yes", ok="no", ok="no")
AC_MSG_RESULT($ok)
if test "$ok" = "no"; then
	AC_DEFINE(NO_FORTRAN_FUNCTIONS,1,[Define if calling Fortran functions directly doesn't work.])
fi

##############################################################################

# grab the various EIGS_ flags from eigensolver.h and substitute them
# into mpb.scm.in (are not regular expressions fun?):
# (note that autoconf gets unhappy if we have newlines; we must delete them)
EIGS_FLAGS_SCM=`grep '#define EIGS_' $srcdir/src/matrices/eigensolver.h |sed 's/#define \([A-Za-z0-9_]*\) (1<<\([0-9]*\))/(define \1 (pow2 \2))/' |sed 's/#define \([A-Za-z0-9_]*\) (\(.*\))/(define \1 (+ \2))/' |sed 's/ *| */ /g' | tr -d '\n'`
AC_SUBST(EIGS_FLAGS_SCM)

##############################################################################

# Find a shell that supports job control, for mpb-split:
AC_PATH_PROG(SH_JOB, bash ksh ash, sh)

if test x"$SH_JOB" = xsh; then
	AC_MSG_WARN([Couldn't find bash or ksh shell on your system.  The mpb-split command may not work if your sh doesn't support job control.])
fi

##############################################################################
# The following function is used only for debugging.  Note that
# we must test for it *after* setting the compiler flags (which
# affect whether it is declared, as it is a GNU extension).
# We need to #include <stdio.h> because that somehow affects whether
# the function is declared with older gcc versions.  We need
# to use AC_TRY_COMPILE because the test in AC_HAVE_DECL seems
# to be optimized out.

AC_CHECK_FUNCS(feenableexcept) dnl GNU libc fp exception control function
AC_MSG_CHECKING([whether feenableexcept declaration is usable])
feenableexcept_decl_ok=yes
AC_TRY_COMPILE([#include <stdio.h>
#define _GNU_SOURCE 1
#include <fenv.h>], [feenableexcept(0);],
[AC_DEFINE([HAVE_DECL_FEENABLEEXCEPT],[1],[Define if fenv.h declares this.])],
[feenableexcept_decl_ok=no])
AC_MSG_RESULT($feenableexcept_decl_ok)

##############################################################################
# The following function is used only for debugging.  Note that
# we must test for it *after* setting the compiler flags (which
# affect whether it is declared, as it is a GNU extension).
# We need to #include <stdio.h> because that somehow affects whether
# the function is declared with older gcc versions.  We need
# to use AC_TRY_COMPILE because the test in AC_HAVE_DECL seems
# to be optimized out.

AC_CHECK_FUNCS(feenableexcept) dnl GNU libc fp exception control function
AC_MSG_CHECKING([whether feenableexcept declaration is usable])
feenableexcept_decl_ok=yes
AC_TRY_COMPILE([#include <stdio.h>
#define _GNU_SOURCE 1
#include <fenv.h>], [feenableexcept(0);],
[AC_DEFINE([HAVE_DECL_FEENABLEEXCEPT],[1],[Define if fenv.h declares this.])],
[feenableexcept_decl_ok=no])
AC_MSG_RESULT($feenableexcept_decl_ok)

##############################################################################
AC_SUBST(MPB_SUFFIX)
program_transform_name="s,_,-,g;$program_transform_name"
##############################################################################

AC_CONFIG_FILES([
	Makefile
	examples/Makefile
	tests/Makefile
	src/Makefile
	src/util/Makefile
	src/matrices/Makefile
	src/matrixio/Makefile
	src/maxwell/Makefile
	mpb/Makefile
	mpb/mpb.scm
	mpb/mpb-split-preinstall
	utils/Makefile
	meb/Makefile
	meb/meb.scm
	meb/meb-split-preinstall
])
AC_OUTPUT
