---
# Installation
---

**Note**: Installing MPB from source can be challenging for novice users. As a simple workaround, the latest version of MPB preinstalled on Ubuntu can be accessed on Amazon Web Services (AWS) Elastic Compute Cloud (EC2) as a free [Amazon Machine Image (AMI)](https://aws.amazon.com/marketplace/pp/B01KHWH0AS). To access this AMI, follow these [instructions](http://www.simpetuscloud.com/launchsims.html).

In this section, we outline the procedure for installing MPB. Mainly, this consists of downloading and installing various prerequisites. As much as possible, we have attempted to take advantage of existing packages such as BLAS, LAPACK, FFTW, and Guile, in order to make our code smaller, more robust, faster, and more flexible. Unfortunately, this may make the installation of MPB more complicated if you do not already have these packages.

You will also need an ANSI C compiler (gcc is fine) and installation will be easiest on a UNIX-like system (Linux is fine). In the following list, some of the packages are dependent upon packages listed earlier, so you should install them in more-or-less the order given.

Many of these libraries may be available in precompiled binary form, especially for Linux systems. Be aware, however, that library binary packages often come in two parts, `library` and `library-dev`, and *both* are required to compile programs using it.

It is important that you use the *same Fortran compiler* to compile Fortran libraries (like LAPACK) and for configuring MPB. Different Fortran compilers often have incompatible linking schemes. The Fortran compiler for MPB can be set via the `F77` environment variable.

[TOC]

Installation on macOS
---------------------

See the [installation guide for Meep on macOS](https://meep.readthedocs.io/en/latest/Installation/#installation-on-macos) for the easiest way to do this.

Unix Installation Basics
------------------------

### Installation Paths

First, let's review some important information about installing software on Unix systems, especially in regards to installing software in non-standard locations. None of these issues are specific to MPB, but they've caused a lot of confusion among users.

Most of the software below, including MPB, installs under `/usr/local` by default. That is, libraries go in `/usr/local/lib`, programs in `/usr/local/bin`, etc. If you don't have `root` privileges on your machine, you may need to install somewhere else, e.g. under `$HOME/install` (the `install/` subdirectory of your home directory). Most of the programs below use a GNU-style `configure` script, which means that all you would do to install there would be:

```
 ./configure --prefix=$HOME/install
```

when configuring the program. The directories `$HOME/install/lib` etc. are created automatically as needed.

#### Paths for Configuring

There are two further complications. First, if you install in a non-standard location and `/usr/local` is considered non-standard by some proprietary compilers, you will need to tell the compilers where to find the libraries and header files that you installed. You do this by setting two environment variables:

```
 setenv LDFLAGS "-L/usr/local/lib"
 setenv CPPFLAGS "-I/usr/local/include"
```

Of course, substitute whatever installation directory you used. Do this **before** you run the `configure` scripts, etcetera. You may need to include multiple `-L` and `-I` flags separated by spaces if your machine has stuff installed in several non-standard locations. Bourne shell users (e.g. `bash` or `ksh`) should use the `export FOO=bar` syntax instead of `csh`'s `setenv FOO bar`, of course.

You might also need to update your `PATH` so that you can run the executables you installed although `/usr/local/bin/` is in the default `PATH` on many systems. e.g. if we installed in our home directory as described above, we would do:

```
 setenv PATH "$HOME/install/bin:$PATH"
```

#### Paths for Running (Shared Libraries)

Second, many of the packages installed below (e.g. Guile) are installed as shared libraries. You need to make sure that your runtime linker knows where to find these shared libraries. The bad news is that every operating system does this in a slightly different way. The good news is that, when you run `make install` for the packages involving shared libraries, the output includes the necessary instructions specific to your system, so pay close attention! It will say something like `add LIBDIR to the <foobar> environment variable`, where `LIBDIR` will be your library installation directory (e.g. `/usr/local/lib`) and `<foobar>` is some environment variable specific to your system (e.g. `LD_LIBRARY_PATH` on some systems, including Linux). For example, you might do:

```
 setenv LD_LIBRARY_PATH "/usr/local/lib:$LD_LIBRARY_PATH"
```

Note that we just add to the library path variable, and don't replace it in case it contains stuff already. If you use Linux and have `root` privileges, you can instead simply run `/sbin/ldconfig`, first making sure that a line `/usr/local/lib` (or whatever) is in `/etc/ld.so.conf`.

If you don't want to type these commands every time you log in, you can put them in your `~/.cshrc` file or `~/.profile`, or `~/.bash_profile`, depending on your shell.

### Fun with Fortran

MPB, along with many of the libraries it calls, is written in C or C++, but it also calls libraries such as BLAS and LAPACK (see below) that are usually compiled from Fortran. This can cause some added difficulty because of the various linking schemes used by Fortran compilers. Our `configure` script attempts to detect the Fortran linking scheme automatically, but in order for this to work ''you must use the same Fortran compiler and options with MPB as were used to compile BLAS/LAPACK''.

By default, MPB looks for a vendor Fortran compiler first (`f77`, `xlf`, etcetera) and then looks for GNU `g77`. In order to manually specify a Fortran compiler `foobar` you would configure it with `./configure F77=foobar ...`.

If, when you compiled BLAS/LAPACK, you used compiler options that alter the linking scheme (e.g. `g77`'s `-fcase-upper` or `-fno-underscoring`), you will need to pass the same flags to MPB via `./configure FFLAGS=...flags... ...`.

### Picking a Compiler

It is often important to be consistent about which compiler you employ. This is especially true for C++ software. To specify a particular C compiler `foo`, configure with `./configure CC=foo`; to specify a particular C++ compiler `foo++`, configure with `./configure CXX=foo++`; to specify a particular Fortran compiler `foo90`, configure with `./configure F77=foo90`.

### Linux and BSD Binary Packages

If you are installing on your personal Linux or BSD machine, then precompiled binary packages are likely to be available for many of these packages, and may even have been included with your system. On Debian systems, the packages are in `.deb` format and the built-in `apt-get` program can fetch them from a central repository. On Red Hat, SuSE, and most other Linux-based systems, binary packages are in RPM format.  OpenBSD has its "ports" system, and so on.

**Do not compile something from source if an official binary package is available.**  For one thing, you're just creating pain for yourself.  Worse, the binary package may already be installed, in which case installing a different version from source will just cause trouble.

One thing to watch out for is that libraries like LAPACK, Guile, HDF5, etcetera, will often come split up into two or more packages: e.g. a `guile` package and a `guile-devel` package. You need to install **both** of these to compile software using the library.

To build the latest version of MPB from source on Ubuntu 16.04, follow these [instructions](http://www.mail-archive.com/mpb-discuss@ab-initio.mit.edu/msg01039.html).


BLAS and LAPACK
---------------

MPB requires the BLAS and LAPACK libraries for matrix computations.

### BLAS

The first thing you must have on your system is a BLAS implementation. "BLAS" stands for "Basic Linear Algebra Subroutines," and is a standard interface for operations like matrix multiplication. It is designed as a building-block for other linear-algebra applications, and is used both directly by our code and in LAPACK (see below). By using it, we can take advantage of many highly-optimized implementations of these operations that have been written to the BLAS interface. Note that you will need implementations of BLAS levels 1-3.

You can find more BLAS information, as well as a basic implementation, on the [BLAS Homepage](http://www.netlib.org/blas/). Once you get things working with the basic BLAS implementation, it might be a good idea to try and find a more optimized BLAS code for your hardware. Vendor-optimized BLAS implementations are available as part of the Intel MKL, HP CXML, IBM ESSL, SGI sgimath, and other libraries. An excellent, high-performance, free-software BLAS implementation is  [OpenBLAS](http://www.openblas.net). Another is [ATLAS](http://math-atlas.sourceforge.net/).

Note that the generic BLAS does not come with a `Makefile`; compile it with something like: </nowiki>

```
  wget http://www.netlib.org/blas/blas.tgz
  gunzip blas.tgz
  tar xf blas.tar
  cd BLAS
  f77 -c -O3 *.f   # compile all of the .f files to produce .o files
  ar rv libblas.a *.o    #  combine the .o files into a library
  su -c "cp libblas.a /usr/local/lib"   # switch to root and install
```

Replace `-O3` with your favorite optimization options. On Linux, this could be `g77 -O3 -fomit-frame-pointer -funroll-loops -malign-double`. Note that MPB looks for the standard BLAS library with `-lblas`, so the library file should be called `libblas.a` and reside in a standard directory like `/usr/local/lib`. See also below for the `--with-blas=lib` option to MPB's `configure` script, to manually specify a library location.

### LAPACK

LAPACK, the Linear Algebra PACKage, is a standard collection of routines, built on BLAS, for more-complicated (dense) linear algebra operations like matrix inversion and diagonalization. You can download LAPACK from the [LAPACK Home Page](http://www.netlib.org/lapack).

Note that MPB looks for LAPACK by linking with `-llapack`. This means that the library must be called `liblapack.a` and be installed in a standard directory like `/usr/local/lib`. Alternatively, you can specify another directory via the `LDFLAGS` environment variable as described earlier. See also below for the `--with-lapack=''lib''` option to our `configure` script, to manually specify a library location.

We currently recommend installing OpenBLAS which includes LAPACK so you do not need to install it separately.

MPI (parallel machines)
------------------------

Optionally, MPB is able to run on a distributed-memory parallel machine, and to do this we use the standard message-passing interface (MPI). You can learn about MPI from its [homepage](http://www-unix.mcs.anl.gov/mpi/). Most commercial supercomputers already have an MPI implementation installed. The recommended implementation is [Open MPI](http://www.open-mpi.org/). MPI is **not required** to compile the serial version of MPB.

In order for the MPI version to run successfully, we have a slightly nonstandard requirement: each process must be able to read from the disk. This way, Guile can boot for each process and they can all read your control file in parallel. Most commercial supercomputers satisfy this requirement.

If you use MPB with MPI, you should compile HDF5 with MPI support as well. See below.

Also, in order to get good performance, you'll need fast interconnect hardware such as Gigabit Ethernet, InfiniBand, or Myrinet. The speed bottleneck comes from the FFT, which requires most of the communications in MPB. FFTW's MPI transforms ([see below](Installation.md#fftw)) come with benchmark programs that will give you a good idea of whether you can get speedups on your system. Of course, even with slow communications, you can still benefit from the memory savings per CPU for large problems.

As described below, when you configure MPB with MPI support (`--with-mpi`), it installs itself as `mpb-mpi`. See also the [User Interface](Scheme_User_Interface.md) for information on using MPB on parallel machines. Normally, you should *also* install the serial version of MPB, if only to get the `mpb-data` utility, which is not installed with the MPI version.

HDF5 (recommended)
-----------------------------------

We require a portable, standard binary format for outputting the electromagnetic fields and similar volumetric data, and for this we use HDF. If you don't have HDF5, you can still compile MPB, but you won't be able to output the fields or the dielectric function.

HDF is a widely-used, free, portable library and file format for multi-dimensional scientific data, developed in the National Center for Supercomputing Applications (NCSA) at the University of Illinois. You can get HDF and learn about it on the [HDF Home Page](http://www.hdfgroup.org).

We require HDF5 which is supported by a number scientific of visualization tools including our own [h5utils](https://github.com/stevengj/h5utils) utilities.

HDF5 includes parallel I/O support under MPI, which can be enabled by configuring it with `--enable-parallel`. You may also have to set the `CC` environment variable to `mpicc`. Unfortunately, the parallel HDF5 library then does not work with serial code, so you have may have to choose one or the other.

We have some hacks in MPB so that it can do parallel I/O even with the serial HDF5 library. These hacks work okay when you are using a small number of processors, but on large supercomputers we strongly recommend using the parallel HDF5.

**Note:** If you have a version of HDF5 compiled with MPI parallel I/O support, then you need to use the MPI compilers to link to it, even when you are compiling the serial version of MPB.  Just use `./configure CC=mpicc CXX=mpic++` or whatever your MPI compilers are when configuring.

FFTW
----

FFTW is a self-optimizing, portable, high-performance FFT implementation, including both serial and parallel FFTs. You can download FFTW and find out more about it from the [FFTW Home Page](http://www.fftw.org).

If you want to use MPB on a parallel machine with MPI, you will also need to install the MPI FFTW libraries. This just means including `--enable-mpi` in the FFTW `configure` flags.

Readline (optional)
-------------------------

Readline is a library to provide command-line history, tab-completion, emacs keybindings, and other shell-like niceties to command-line programs. This is an *optional* package, but one that can be used by Guile (see below) if it is installed. We recommend installing it. You can download Readline from its [ftp site](ftp://ftp.gnu.org/gnu/readline). Readline is typically preinstalled on Linux systems.

Guile
-----------------------

Guile is required in order to use the Scheme interface, and is strongly recommended. If you don't install it, you can only use the C++ interface.

Guile is an extension/scripting language implementation based on Scheme, and we use it to provide a rich, fully-programmable user interface with minimal effort. It's free, of course, and you can download it from the [Guile Home Page](http://www.gnu.org/software/guile/). Guile is typically included with Linux systems.

- **Important:** Most Linux distributions come with Guile already installed. You can check by seeing whether you can run `guile --version` from the command line. In that case, do **not** install your own version of Guile from source &mdash; having two versions of Guile on the same system will cause problems. However, by default most distributions install only the Guile libraries and not the programming headers &mdash; to compile libctl and MPB, you should install the **guile-devel** or **guile-dev** package.

Autoconf (optional)
-------------------------

If you want to be a developer of the MPB package as opposed to merely a user, you will also need the Autoconf program. Autoconf is a portability tool that generates `configure` scripts to automatically detect the capabilities of a system and configure a package accordingly. You can find out more at the [Autoconf Home Page](https://www.gnu.org/software/autoconf/autoconf.html). `autoconf` is typically installed by default on Linux systems. In order to install Autoconf, you will also need the GNU `m4` program if you do not already have it. See the [m4 Home Page](https://www.gnu.org/software/m4/m4.html).

libctl
------

[libctl](https://libctl.readthedocs.io), which requires Guile, is required to use the Scheme interface, and is strongly r
ecommended. If you don't install it, you can only use the C++ interface. libctl version **3.2 or later** is required.

Instead of using Guile directly, we separated much of the user interface code into a package called libctl, in the hope that this might be
 more generally useful. libctl automatically handles the communication between the program and Guile, converting complicated data structur
 es and so on, to make it even easier to use Guile to control scientific applications. Download libctl from the [libctl page](https://libctl.readthedocs.io), unpack it, and run the usual `configure`, `make`, `make install` sequence. You'll also want to browse the [libctl manual](https://libctl.readthedocs.io), as this will give you a general overview of what the user interface will be like.

If you are not the system administrator of your machine, and/or want to install libctl somewhere else like your home directory, you can do so with the standard `--prefix=dir` option to `configure`. The default prefix is `/usr/local`. In this case, however, you'll need to specify the location of the libctl shared files for MPB, using the `--with-libctl=dir/share/libctl` option to our `configure` script.

MPB
------------------

If you've made it all the way here, you're ready to install the MPB package and start cranking out eigenmodes. You can obtain the latest version from the [Download page](Download.md). Once you've unpacked it, just run:

```
./configure
make
```

to configure and compile the package. See below to install. Hopefully, the `configure` script will correctly detect the BLAS, FFTW, etcetera libraries which have been installed, as well as the C compiler and so on, and the `make` compilation will proceed without a hitch. If not, `configure` accepts several flags to help control its behavior. Some of these are standard, like `--prefix=`*`dir`* to specify an installation directory prefix, and some of them are specific to the MPB package. Use `./configure --help` for more info. The `configure` flags specific to MPB are:

**`--with-inv-symmetry`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Assume [inversion symmetry](Scheme_User_Interface.md#inversion-symmetry) in the dielectric function, allowing us to use real fields in Fourier space instead of complex fields. This gives a factor of 2 benefit in speed and memory. In this case, the MPB program will be installed as `mpbi` instead of `mpb`, so that you can have versions both with and without inversion symmetry installed at the same time. To install *both* `mpb` and `mpbi`, you should do:

```
./configure
make
sudo make install
make distclean
./configure --with-inv-symmetry
make
sudo make install
```

**`--with-hermitian-eps`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Support the use of [complex-hermitian dielectric tensors](Scheme_User_Interface.md#material-type) corresponding to magnetic materials, which break inversion symmetry.

**`--enable-single`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Use single precision (C `float`) instead of the default double precision (C `double`) for computations. Not recommended.

**`--without-hdf5`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Don't use the HDF5 library for field and dielectric function output. In which case, no field output is possible.

**`--with-mpi`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Attempt to compile a parallel version of MPB using MPI; the resulting program will be installed as `mpb-mpi`. Requires [MPI](Installation.md#mpi-parallel-machines) and [MPI FFTW](Installation.md#fftw) libraries to be installed, as described above.

Does *not* compile the serial MPB, or `mpb-data`; if you want those, you should `make distclean` and compile/install them separately.

**`--with-mpi`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Can be used along with `--with-inv-symmetry`, in which case the program is installed as `mpbi-mpi`.

**`--with-openmp`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Attempt to compile a shared-memory parallel version of MPB using OpenMP. The resulting program will be installed as `mpb` and FFTs will use OpenMP parallelism. Requires OpenMP FFTW libraries to be installed.

**`--with-libctl=dir`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
If libctl was installed in a nonstandard location (i.e. neither `/usr` nor `/usr/local`), you need to specify the location of the libctl directory, `dir`. This is either `prefix/share/libctl`, where `prefix` is the installation prefix of libctl, or the original libctl source code directory.

**`--with-blas=lib`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The `configure` script automatically attempts to detect accelerated BLAS libraries, like DXML (DEC/Alpha), SCSL and SGIMATH (SGI/MIPS), ESSL (IBM/PowerPC), ATLAS, and PHiPACK. You can, however, force a specific library name to try via `--with-blas=lib`.

**`--with-lapack=lib`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Cause the `configure` script to look for a LAPACK library called `lib`. The default is to use `-llapack`.

**`--disable-checks`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Disable runtime checks. Not recommended. The disabled checks shouldn't take up a significant amount of time anyway.

**`--enable-prof`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Compile for performance profiling.

**`--enable-debug`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Compile for debugging, adding extra runtime checks and so on.

**`--enable-debug-malloc`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Use special memory-allocation routines for extra debugging (to check for array overwrites, memory leaks, etcetera).

**`--with-efence`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
More debugging: use the [Electric Fence](http://elinux.org/Electric_Fence) library, if available, for extra runtime array bounds-checking.

You can further control `configure` by setting various environment variables, such as:

-   `CC`: the C compiler command
-   `CFLAGS`: the C compiler flags (defaults to `-O3`).
-   `CPPFLAGS`: `-I`*`dir`* flags to tell the C compiler additional places to look for header files.
-   `LDFLAGS`: `-L`*`dir`* flags to tell the linker additional places to look for libraries.
-   `LIBS`: additional libraries to link against.

Once compiled, the main program, as opposed to various test programs, resides in the `mpb-ctl/` subdirectory, and is called `mpb`. You can install this program under `/usr/local` or elsewhere, if you used the `--prefix` flag for `configure`, by running:

```
sudo make install
```

The "sudo" command is to switch to `root` for installation into system directories. You can just do `make install` if you are installing into your home directory instead.

If you make a mistake (e.g. you forget to specify a needed `-L`*` dir`* flag) or in general want to start over from a clean slate, you can restore MPB to a pristine state by running:

```
make distclean
```
