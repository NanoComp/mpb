---
Developer Information
---

Here, we begin with a brief overview of what the program is computing, and then describe how the program and computation are broken up into different portions of the code.

A [Chinese version](http://ab-initio.mit.edu/mpb/mpb-devel-zh_CN.pdf) of this document is also available.

[TOC]

The Mathematics of MPB
----------------------

This section provides a whirlwind tour of the mathematics of photonic band structure calculations and the algorithms that we employ. For more detailed information, see:

-   [Photonic Crystals: Molding the Flow of Light](http://ab-initio.mit.edu/book), by J. D. Joannopoulos, S. G. Johnson, R. D. Meade, and J. N. Winn (Princeton, 2008).
-   Steven G. Johnson and J. D. Joannopoulos, [Block-iterative frequency-domain methods for Maxwell's equations in a planewave basis](http://www.opticsexpress.org/abstract.cfm?URI=OPEX-8-3-173), *Optics Express* **8**, no. 3, 173-190 (2001).

MPB takes a periodic dielectric structure and computes the *eigenmodes* of that structure, which are the electromagnetic waves that can propagate through the structure with a definite frequency. This corresponds to solving an eigenvalue problem

\[\hat\Theta \mathbf{H} = \frac{\omega^2}{c^2} \mathbf{H},\] where \(\mathbf{H}\) is the magnetic field, \(\omega\) is the frequency, and \(\hat\Theta\) is the Maxwell operator

\[\hat\Theta = \nabla\times \frac{1}{\varepsilon} \nabla\times \,.\] We also have an additional constraint, that \(\nabla \cdot \mathbf{H}\) be zero (the magnetic field must be "transverse").

Since the structure is periodic, we can also invoke Bloch's theorem to write the states in the form:

\[\mathbf{H} = \mathbf{H}_\mathbf{k}(\mathbf{x}) e^{i \mathbf{k} \cdot \mathbf{x}},\] where \(\mathbf{H}_\mathbf{k}\) is a periodic function (the Bloch envelope) and \(\mathbf{k}\) is the Bloch wavevector. So, at each k-point (Bloch wavevector), we need to solve for a discrete set of eigenstates, the photonic bands of the structure.

To solve for the eigenstates on a computer, we must expand the magnetic field in some basis, where we truncate the basis to some finite number of points to discretize the problem. For example, we could use a traditional finite-element basis in which the field is taken on a finite number of mesh points and linearly interpolated in between. However, it is expensive to enforce the transversality constraint in this basis. Instead, we use a Fourier (spectral) basis, expanding the periodic part of the field as a sum of planewaves:

\[\mathbf{H}_\mathbf{k} = \sum_\mathbf{G} \mathbf{h}_\mathbf{G} e^{i \mathbf{G} \cdot \mathbf{x}} \, .\] In this basis, the transversality constraint is easy to maintain, as it merely implies that the planewave amplitudes \(\mathbf{h}_\mathbf{G}\) must be orthogonal to \(\mathbf{k}+\mathbf{G}\).

In order to find the eigenfunctions, we could compute the elements of \(\hat\Theta\) explicitly in our basis, and then call LAPACK or some similar code to find the eigenvectors and eigenvalues. For a three-dimensional calculation, this could mean finding the eigenvectors of a matrix with millions of elements on a side--daunting merely to store, much less compute. Fortunately, we only want to know a few eigenvectors, not hundreds of thousands, so we can use much less expensive *iterative* methods that don't require us to store \(\hat\Theta\) explicitly.

Iterative eigensolvers require only that one supply a routine to operate \(\hat\Theta\) on a vector (function). Starting with an initial guess for the eigenvector, they then converge quickly to the actual eigenvector, stopping when the desired tolerance is achieved. There are many iterative eigensolver methods; we use a preconditioned block minimization of the Rayleigh quotient which is further described in the file `src/matrices/eigensolver.c`. In the Fourier basis, applying \(\hat\Theta\) to a function is relatively easy: the curls become cross products with \(i(\mathbf{k}+\mathbf{G})\); the multiplication by \(1/\varepsilon\) is performed by using an [FFT](https://en.wikipedia.org/wiki/Fast_Fourier_transform) to transform to the spatial domain, multiplying, and then transforming back with an inverse FFT. For more information and references on iterative eigensolvers, see the paper cited above.

We also support a "targeted" eigensolver. A typical iterative eigensolver finds the *p* lowest eigenvalues and eigenvectors. Instead, we can find the *p* eigenvalues closest to a given frequency \(\omega_0\) by solving for the eigenvalues of \((\hat\Theta-\omega_0^2/c^2)^2\) instead of \(\hat\Theta\). This new operator has the same eigenvectors as \(\hat\Theta\), but its eigenvalues have been shifted to make those closest to \(\omega_0\) the smallest. This is not really the best algorithm to find interior eigenvalues like this; a future version of MPB may use ARPACK-style shift-and-invert Arnoldi, or perhaps the Jacobi-Davidson algorithm.

The eigensolver we use is preconditioned, which means that convergence can be greatly improved by suppling a good preconditioner matrix. Finding a good preconditioner involves making an approximate inverse of \(\hat\Theta\), and is something of a black art with lots of trial and error.

Dielectric Function Computation
-------------------------------

The initialization of the dielectric function deserves some additional discussion, both because it is crucial for good convergence, and because we use somewhat complicated algorithms for performance reasons.

To ameliorate the convergence problems caused in a planewave basis by a discontinuous dielectric function, the dielectric function is smoothed (averaged) at the resolution of the grid. Another way of thinking about it is that this brings the average dielectric constant (over the grid) closer to its true value. Since different polarizations of the field prefer different averaging methods, one has to construct an effective dielectric tensor at the boundaries between dielectrics, as described by the paper referenced above.

This averaging has two components. First, at each grid point the dielectric constant ($\varepsilon$) and its inverse are averaged over a uniform mesh extending halfway to the neighboring grid points. The mesh resolution is controlled by the `mesh-size` user input variable. Second, for grid points on the boundary between two dielectrics, we compute the vector normal to the dielectric interface; this is done by averaging the "dipole moment" of the dielectric function over a spherically-symmetric distribution of points. The normal vector and the two averages of epsilon are then combined into an effective dielectric tensor for the grid point.

All of this averaging is handled by a subroutine in `src/maxwell/` (see below) that takes as input a function $\varepsilon$(**r**), which returns the dielectric constant for a given position **r**. This epsilon function must be as efficient as possible, because it is evaluated a large number of times: the size of the grid multiplied by `mesh-size`<sup>3</sup> (in three dimensions).

To specify the geometry, the user provides a list of geometric objects (blocks, spheres, cylinders and so on). These are parsed into an efficient data structure and are used to to provide the epsilon function described above. All of this is handled by the libctlgeom component of libctl, described below. At the heart of the epsilon function is a routine to return the geometric object enclosing a given point, taking into account the fact that the objects are periodic in the lattice vectors. Our first algorithm for doing this was a simple linear search through the list of objects and their translations by the lattice vectors, but this proved to be too slow, especially in supercell calculations where there are many objects. We addressed the performance problem in two ways. First, for each object we construct a bounding box, with which point inclusion can be tested rapidly. Second, we build a hierarchical tree of bounding boxes, recursively partitioning the set of objects in the cell. This allows us to search for the object containing a point in a time logarithmic in the number of objects instead of linear as before.

Code Organization
-----------------

The code is organized to keep the core computation independent of the user interface, and to keep the eigensolver routines independent of the operator they are computing the eigenvector of. The computational code is located in the `src/` directory, with a few major subdirectories, described below. The Guile-based user interface is completely contained within the `mpb-ctl/` directory.

### src/matrices/

This directory contains the eigensolver, in `eigensolver.c`, to which you pass an operator and it returns the eigenvectors. Eigenvectors are stored using the `evectmatrix` data structure, which holds `p` eigenvectors of length `n`, potentially distributed over `n` in MPI. See `src/matrices/README` for more information about the data structures. In particular, you should use the supplied functions (`create_evectmatrix`, etcetera) to create and manipulate the data structures, where possible.

The type of the eigenvector elements is determined by `scalar.h`, which sets whether they are real or complex and single or double precision. This is, in turn, controlled by the `--disable-complex` and `--enable-single` parameters to the `configure` script at install-time. `scalar.h` contains macros to make it easier to support both real and complex numbers elsewhere in the code.

Also in this directory is `blasglue.c`, a set of wrapper routines to make it convienient to call BLAS and LAPACK routines from C instead of Fortran.

### src/util/

As its name implies, this is simply a number of utility routines for use elsewhere in the code. Of particular note is `check.h`, which defines a `CHECK(condition, error-message)` macro that is used extensively in the code to improve robustness. There are also debugging versions of malloc/free (which perform lots of paranoia tests, enabled by `--enable-debug-malloc` in `configure`), and MPI glue routines that allow the program to operate without the MPI libraries.

### src/matrixio

This section contains code to abstract I/O for eigenvectors and similar matrices, providing a simpler layer on top of the HDF5 interface. This could be modified to support other I/O formats.

### src/maxwell/

The `maxwell/` directory contains all knowledge of Maxwell's equations used by the program. It implements functions to apply the Maxwell operator to a vector (in `maxwell_op.c`) and compute a good preconditioner (in `maxwell_pre.c`). These functions operate upon a representation of the fields in a transverse Fourier basis.

In order to use these functions, one must first initialize a `maxwell_data` structure with `create_maxwell_data` (defined in `maxwell.c`) and specify a k point with `update_maxwell_data_k`. One must also initialize the dielectric function using `set_maxwell_dielectric` by supplying a function that returns the dielectric constant for any given coordinate. You can also restrict yourself to TE or TM polarizations in two dimensions by calling `set_maxwell_data_polarization`.

This directory also contains functions `maxwell_compute_dfield`, etcetera, to compute the position-space fields from the Fourier-transform representation returned by the eigensolver.

### mpb-ctl/

Here is the Guile-based user interface code for the eigensolver. Instead of using Guile directly, this code is built on top of the `libctl` library as described in previous sections. This means that the user-interface code (in `mpb.c`) is fairly short, consisting of a number of small functions that are callable by the user from Guile.

The core of the user interface is the file `mpb.scm`, the *specifications file* for libctl as described in the [libctl manual](http://ab-initio.mit.edu/libctl/doc/). Actually, `mpb.scm` is generated by `configure` from `mpb.scm.in` in order to substitute in parameters like the location of the libctl library. You should only edit `mpb.scm.in` directly. You can regenerate `mpb.scm` simply by running `./config.status` instead of re-running `configure`.

The specifications file defines the data structures and subroutines that are visible to the Guile user. It also defines a number of Scheme subroutines for the user to call directly, like `(run)`. It is often simpler and more flexible to define functions like this in Scheme rather than in C.

All of the code to handle the geometric objects resides in libctlgeom, a set of Scheme and C utility functions included with libctl (see the file `utils/README` in the libctl package). These functions could also be useful in other programs, such as a time-domain Maxwell's equation simulator.