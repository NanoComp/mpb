# Process this file with autoconf to produce a configure script.
AC_INIT
AC_CONFIG_SRCDIR([src/matrices/matrices.c])
AC_CONFIG_HEADER(src/config.h)

MPB_VERSION=1.3
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

# Checks for programs.
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_PROG_MAKE_SET

# check how to transform the name of the installed program:
AC_ARG_PROGRAM

AC_ARG_ENABLE(single, [AC_HELP_STRING([--enable-single],[compile for single precision])], enable_single=$enableval, enable_single=no)
if test "$enable_single" = "yes"; then
        AC_DEFINE(SCALAR_SINGLE_PREC,1,[Define to use single precision])
	echo "*********************** SINGLE PRECISION ***********************"
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

BUILD_MPB_DATA=mpb-data  # mpb-data program to build, if any
AC_SUBST(BUILD_MPB_DATA)

INSTALL_MPB_SPLIT=install-mpb-split  # mpb-split target to install, if any
AC_SUBST(INSTALL_MPB_SPLIT)

AC_ARG_WITH(inv-symmetry,[AC_HELP_STRING([--with-inv-symmetry],[take advantage of (and require) inv. sym.])], inv_sym=$withval, inv_sym=no)
if test "$inv_sym" = "no"; then
        AC_DEFINE(SCALAR_COMPLEX,1,[Define to use complex fields and
	                            not to require inversion symmetry])
	SCALAR_TYPE=complex
else
	SCALAR_TYPE=real
	# Make installed program mpbi instead of mpb, but don't change mpb-data
	program_transform_name='s,^mpb$$,mpbi,'";$program_transform_name"
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
# Check for FFTW libraries
if test "$enable_single" = "yes"; then
	AC_CHECK_LIB(sfftw, fftw)
	if test x"$inv_sym" = xyes; then
		AC_CHECK_LIB(srfftw, rfftw)
	fi
else
	AC_CHECK_LIB(dfftw, fftw)
	if test x"$inv_sym" = xyes; then
		AC_CHECK_LIB(drfftw, rfftw)
	fi
fi

if test x = x"`echo $LIBS | egrep 'l[[ds]]fftw'`"; then
	AC_CHECK_LIB(fftw, fftw)
	if test x"$inv_sym" = xyes; then
		AC_CHECK_LIB(rfftw, rfftw)
	fi
fi

if test x = x"`echo $LIBS | egrep 'l[[ds]]*fftw'`"; then
	AC_MSG_ERROR([The FFTW libraries could not be found.  Make sure FFTW is installed, and that LDFLAGS and CPPFLAGS are set appropriately if necessary.])
fi

if test x"$inv_sym" = xyes; then
	if test x = x"`echo $LIBS | egrep 'l[[ds]]*rfftw'`"; then
		AC_MSG_ERROR([The RFFTW libraries, which are required to compile MPB with inversion symmetry, could not be found.  These should have been installed as a part of FFTW.])
	fi
fi

##############################################################################
# Checks for BLAS and LAPACK libraries:

sinclude(acx_blas.m4)
sinclude(acx_lapack.m4)

ACX_BLAS([], [AC_MSG_ERROR([BLAS was not found!])])
ACX_LAPACK([], [AC_MSG_ERROR([LAPACK was not found!])])

LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

##############################################################################
# Check for MPI library

AC_ARG_WITH(mpi, [AC_HELP_STRING([--with-mpi],[enable MPI parallelization])], with_mpi=$withval, with_mpi=no)

if test "x$with_mpi" = "xyes"; then
	AC_CHECK_PROGS(MPICC, mpicc hcc mpcc mpcc_r, $CC)
        CC="$MPICC"
        if test x = x"$MPILIBS"; then
                AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])
        fi
        if test x = x"$MPILIBS"; then
                AC_CHECK_LIB(mpi, MPI_Init, [MPILIBS="-lmpi"])
        fi

        if test x = x"$MPILIBS"; then
                AC_MSG_ERROR(could not find mpi library for --with-mpi)
	fi

	LIBS="$MPILIBS $LIBS"

	# Make installed program mpb_mpi instead of mpb
        if test x"$inv_sym" = xyes; then
	        program_transform_name='s,^mpb$$,mpbi-mpi,'";$program_transform_name"
	else
	        program_transform_name='s,^mpb$$,mpb-mpi,'";$program_transform_name"
	fi

	# Check for fftw_mpi library.
	if test "$enable_single" = "yes"; then
        	AC_CHECK_LIB(sfftw_mpi, fftwnd_mpi)
	        if test x"$inv_sym" = xyes; then
	                AC_CHECK_LIB(srfftw_mpi, rfftwnd_mpi)
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

	BUILD_MPB_DATA="" # don't build/install MPI-linked mpb-data
	INSTALL_MPB_SPLIT="" # don't install mpb-split, since we don't know
	                     # how to run mpb-mpi in general
fi

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
	
	save_LIBS_0="$LIBS" # save, to check later if we found any library
        AC_CHECK_LIB(guile, gh_eval_str)
	if test x"$save_LIBS_0" = x"$LIBS"; then
		AC_MSG_ERROR([Couldn't find the required Guile library.])
	fi
fi

AC_CHECK_FUNCS(scm_make_smob_type)

# Check how to activate Guile readline support:

ACTIVATE_READLINE=""
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

AC_CHECK_HEADER(ctl.h, [], [AC_MSG_ERROR([Couldn't find the libctl header.])])

# Check libctl version >= LIBCTL_MAJOR.LIBCTL_MINOR.LIBCTL_BUGFIX
LIBCTL_MAJOR=2; LIBCTL_MINOR=1; LIBCTL_BUGFIX=1
AC_MSG_CHECKING([whether libctl version is at least ${LIBCTL_MAJOR}.${LIBCTL_MINOR}.${LIBCTL_BUGFIX}])
AC_EGREP_CPP(yes, [[
#include <ctl.h>
#if LIBCTL_MAJOR_VERSION > $LIBCTL_MAJOR || (LIBCTL_MAJOR_VERSION == $LIBCTL_MAJOR && (LIBCTL_MINOR_VERSION > $LIBCTL_MINOR || (LIBCTL_MINOR_VERSION == $LIBCTL_MINOR && LIBCTL_BUGFIX_VERSION >= $LIBCTL_BUGFIX)))
  yes
#endif
]], [AC_MSG_RESULT(ok)], [AC_MSG_ERROR([libctl version ${LIBCTL_MAJOR}.${LIBCTL_MINOR}.${LIBCTL_BUGFIX} or later is required])])

##############################################################################
# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS(unistd.h getopt.h)

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

AC_CONFIG_FILES([Makefile tests/Makefile src/util/Makefile src/matrices/Makefile src/matrixio/Makefile src/maxwell/Makefile mpb-ctl/Makefile mpb-ctl/mpb.scm mpb-ctl/mpb-split utils/Makefile])
AC_OUTPUT
