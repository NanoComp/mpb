---
title: MPB Installation
permalink: /MPB_Installation/
---

In this section, we outline the procedure for installing the MIT Photonic-Bands package. Mainly, this consists of downloading and installing various prerequisites. As much as possible, we have attempted to take advantage of existing packages such as BLAS, LAPACK, FFTW, and GNU Guile, in order to make our code smaller, more robust, faster, and more flexible. Unfortunately, this may make the installation of MPB more complicated if you do not already have these packages.

You will also need an ANSI C compiler, of course (gcc is fine), and installation will be easiest on a UNIX-like system (Linux is fine). In the following list, some of the packages are dependent upon packages listed earlier, so you should install them in more-or-less the order given.

**Note:** Many of these libraries may be available in precompiled binary form, especially for GNU/Linux systems. Be aware, however, that library binary packages often come in two parts, `library` and `library-dev`, and *both* are required to compile programs using it.

**Note:** It is important that you use the *same Fortran compiler* to compile Fortran libraries (like LAPACK) and for configuring MPB. Different Fortran compilers often have incompatible linking schemes. (The Fortran compiler for MPB can be set via the `F77` environment variable.)

**Note:** The latest, pre-installed versions of MPB and Meep running on Ubuntu can also be accessed on Amazon Web Services (AWS) Elastic Compute Cloud (EC2) as a free [Amazon Machine Image (AMI)](https://aws.amazon.com/marketplace/pp/B01KHWH0AS). To access this AMI, follow these [instructions](http://www.simpetuscloud.com/launchsims.html).

Installation on MacOS X
-----------------------

See [Meep installation on OS X](/Meep_Installation#Installation_on_MacOS_X "wikilink") for the easiest way to do this.

BLAS and LAPACK
---------------

MPB requires the BLAS and LAPACK libraries for matrix computations.

MPI *(for parallel MPB)*
------------------------

Optionally, MPB is able to run on a distributed-memory parallel machine, and to do this we use the standard MPI message-passing interface. You can learn about MPI from the [MPI Home Page](http://www-unix.mcs.anl.gov/mpi/). Most commercial supercomputers already have an MPI implementation installed; two free MPI implementations that you can install yourself are [MPICH](http://www-unix.mcs.anl.gov/mpi/mpich/) and [LAM](http://www.lam-mpi.org/); a promising new implementation is [Open MPI](http://www.open-mpi.org/). MPI is *not required* to compile the ordinary, uniprocessor version of our software.

In order for the MPI version of our software to run successfully, we have a slightly nonstandard requirement: each process must be able to read from the disk. (This way, Guile can boot for each process and they can all read your control file in parallel.) Many (most?) commercial supercomputers, Linux [Beowulf](http://www.beowulf.org) clusters, etcetera, satisfy this requirement.

Also, in order to get good performance, you'll need fast interconnect hardware such as Myrinet; 100Mbps Ethernet probably won't cut it. The speed bottleneck comes from the FFT, which requires most of the communications in MPB. FFTW's MPI transforms ([see below](/#FFTW "wikilink")) come with benchmark programs that will give you a good idea of whether you can get speedups on your system. Of course, even with slow communications, you can still benefit from the memory savings per CPU for large problems.

As described [below](/#MIT_Photonic_Bands "wikilink"), when you configure MPB with MPI support (`--with-mpi`), it installs itself as `mpb-mpi`. See also the \[user-ref.html\#mpb-mpi user reference\] section for information on using MPB on parallel machines. Normally, you should *also* install the serial version of MPB, if only to get the `mpb-data` utility, which is not installed with the MPI version.

MPI support in MPB is thanks to generous support from [Clarendon Photonics](http://www.clarendonphotonics.com/).

HDF5 *(optional, strongly advised)*
-----------------------------------

We require a portable, standard binary format for outputting the electromagnetic fields and similar volumetric data, and for this we use HDF. (If you don't have HDF5, you can still compile MPB, but you won't be able to output the fields or the dielectric function.)

FFTW
----

FFTW is a self-optimizing, portable, high-performance FFT implementation, including both serial and parallel FFTs. You can download FFTW and find out more about it from the [FFTW Home Page](http://www.fftw.org).

If you want to use MPB on a parallel machine with MPI, you will also need to install the MPI FFTW libraries (this just means including `--enable-mpi` in the FFTW `configure` flags).

GNU Readline *(optional)*
-------------------------

GNU Readline is a library to provide command-line history, tab-completion, emacs keybindings, and other shell-like niceties to command-line programs. This is an *optional* package, but one that can be used by Guile (see below) if it is installed; we recommend getting it. You can [download Readline from the GNU ftp site](ftp://ftp.gnu.org/gnu/readline). (Readline is typically preinstalled on GNU/Linux systems).

GNU Guile
---------

GNU Autoconf *(optional)*
-------------------------

If you want to be a developer of the MPB package (as opposed to merely a user), you will also need the GNU Autoconf program. Autoconf is a portability tool that generates `configure` scripts to automatically detect the capabilities of a system and configure a package accordingly. You can find out more at the [Autoconf Home Page](http://sources.redhat.com/autoconf/) (autoconf is typically installed by default on Linux systems). In order to install Autoconf, you will also need the GNU `m4` program if you do not already have it (see the [GNU m4 Home Page](http://www.seindal.dk/rene/gnu/)).

libctl
------

MIT Photonic Bands
------------------

Okay, if you've made it all the way here, you're ready to install the MPB package and start cranking out eigenmodes. (You can download the latest version and read this manual at the [MIT Photonic-Bands Homepage](/MIT_Photonic_Bands "wikilink").) Once you've unpacked it, just run:

`./configure`
`make`

to configure and compile the package (see below to install). Hopefully, the `configure` script will correctly detect the BLAS, FFTW, etcetera libraries that you've dutifully installed, as well as the C compiler and so on, and the `make` compilation will proceed without a hitch. If not, it's a [Simple Matter of Programming](http://www.catb.org/jargon/html/S/SMOP.html) to correct the problem. `configure` accepts several flags to help control its behavior. Some of these are standard, like `--prefix=`*`dir`* to specify and installation directory prefix, and some of them are specific to the MPB package (`./configure` `--help` for more info). The `configure` flags specific to MPB are:

`--with-inv-symmetry`
Assume \[user-ref.html\#inv-symmetry inversion symmetry\] in the dielectric function, allowing us to use real fields (in Fourier space) instead of complex fields. This gives a factor of 2 benefit in speed and memory. In this case, the MPB program will be installed as `mpbi` instead of `mpb`, so that you can have versions both with and without inversion symmetry installed at the same time. To install *both* `mpb` and `mpbi`, you should do:

`./configure`
`make`
`su -c "make install"`
`make distclean`
`./configure --with-inv-symmetry`
`make`
`su -c "make install"`

`--with-hermitian-eps`
Support the use of \[user-ref.html\#dielectric-anisotropic complex-hermitian dielectric tensors\] (corresponding to magnetic materials, which break inversion symmetry).

`--enable-single`
Use single precision (C `float`) instead of the default double precision (C `double`) for computations. (Not recommended.)

`--without-hdf5`
Don't use the HDF5 library for field and dielectric function output. (In which case, no field output is possible.)

`--with-mpi`
Attempt to compile a \[user-ref.html\#mpb-mpi parallel version of MPB\] using MPI; the resulting program will be installed as `mpb-mpi`. Requires \[\#mpi MPI\] and \[\#fftw MPI FFTW\] libraries to be installed, as described above.

Does *not* compile the serial MPB, or `mpb-data`; if you want those, you should `make` `distclean` and compile/install them separately.

`--with-mpi` *can* be used along with `--with-inv-symmetry`, in which case the program is installed as `mpbi-mpi` (try typing that five times quickly).

`--with-openmp`
Attempt to compile a shared-memory parallel version of MPB using OpenMP; the resulting program will be installed as `mpb`, and FFTs will use OpenMP parallelism. Requires OpenMP FFTW libraries to be installed.

`--with-libctl=`*`dir`*
If libctl was installed in a nonstandard location (i.e. neither `/usr` nor `/usr/local`), you need to specify the location of the libctl directory, *`dir`*. This is either *`prefix`*`/share/libctl`, where *`prefix`* is the installation prefix of libctl, or the original libctl source code directory.

`--with-blas=`*`lib`*
The `configure` script automatically attempts to detect accelerated BLAS libraries, like DXML (DEC/Alpha), SCSL and SGIMATH (SGI/MIPS), ESSL (IBM/PowerPC), ATLAS, and PHiPACK. You can, however, force a specific library name to try via `--with-blas=`*`lib`*.

`--with-lapack=`*`lib`*
Cause the `configure` script to look for a LAPACK library called *`lib`* (the default is to use `-llapack`).

`--disable-checks`
Disable runtime checks. (Not recommended; the disabled checks shouldn't take up a significant amount of time anyway.)

`--enable-prof`
Compile for performance profiling.

`--enable-debug`
Compile for debugging, adding extra runtime checks and so on.

`--enable-debug-malloc`
Use special memory-allocation routines for extra debugging (to check for array overwrites, memory leaks, etcetera).

`--with-efence`
More debugging: use the [Electric Fence](http://perens.com/FreeSoftware/) library, if available, for extra runtime array bounds-checking.

You can further control `configure` by setting various environment variables, such as:

-   `CC`: the C compiler command
-   `CFLAGS`: the C compiler flags (defaults to `-O3`).
-   `CPPFLAGS`: `-I`*`dir`* flags to tell the C compiler additional places to look for header files.
-   `LDFLAGS`: `-L`*`dir`* flags to tell the linker additional places to look for libraries.
-   `LIBS`: additional libraries to link against.

Once compiled, the main program (as opposed to various test programs) resides in the `mpb-ctl/` subdirectory, and is called `mpb`. You can install this program under /usr/local (or elsewhere, if you used the `--prefix` flag for `configure`), by running:

`su -c "make install"`

The "su" command is to switch to `root` for installation into system directories. You can just do `make` `install` if you are installing into your home directory instead.

If you make a mistake (e.g. you forget to specify a needed `-L`*`dir`* flag) or in general want to start over from a clean slate, you can restore MPB to a pristine state by running:

`make distclean`

[Category:MPB](/Category:MPB "wikilink")