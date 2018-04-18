# MPB Release Notes

## MPB 1.6.2

4/18/2018

  * Fixed handling of 0-argument band functions (like `randomize-fields`)
    in Guile 2.2.

  * Bugfix in `dot-eigenvectors` for μ ≠ 1 (#37).

## MPB 1.6.1

1/18/2018

  * Corrected some minor release glitches.

## MPB 1.6

1/18/2018

  * Support for Guile 2.2.

  * Support for magnetic materials (mu) via `mu`/`mu-diag`/`mu-offdiag`
    material properties (thanks to Ling Lu).

  * libmpb (for use in Meep) can now be built without Guile.

  * `sqmatrix-diagm` and `sqmatrix-eigvals` functions (thanks to Ling Lu).

  * Migrated documentation to github/markdown/readthedocs (#22).

## MPB 1.5

4/2/2014.

  * MPB now also installs a library, for use from within Meep 1.2 or later.

  * Support Guile 2.x.

  * Support FFTW version 3.x (in addition to FFTW 2.x).  Version 3.3 or later
    is required for MPI parallelism.

  * Support for OpenMP parallelism.

  * Use more accurate subpixel averaging algorithm for interfaces between
    anisotropic material (see Kottke et al, PRE 77, 036611, 2008).

  * Use more accurate geometry routines in recent libctl versions (improves
    and speeds subpixel averaging).

  * Support using a different k origin in find-k.

  * Added epsilon-func wrapper for material-func (similar to Meep).

  * Added compute-1-group-velocity to compute group velocity of a single band.

  * Added kinterpolate-uniform function, to interpolate with roughly
    uniform spacing in reciprocal space.

  * Added optimize-grid-size! function to round the grid size to a size
    that can be handled more efficiently.  (Only affects resolution, not
    lattice vectors.)

  * Allow user to set filename-prefix to false to disable its use.

  * Added mpb-data -P <angle> option to change the phase angle of the output.

  * Added compute-field-divergence and routines to get and output the
    bound charge density: get-charge-density and output-charge-density.

  * resolution is now an arbitrary real number, not just an integer,
    although of course MPB eventually computes an integer grid size.

  * Support HDF5 1.8.

  * Fix recurring "non positive-definite matrix in potrf" errors that
    were arising due to roundoff errors preventing matrix inversion.

  * Fix output-at-kpoint to avoid sensitivity to roundoff errors.

  * Bug fix in parallel HDF5 support: HDF5 compiled for parallel I/O (MPI)
    now works.

  * Bug fix in field-map!, thanks to Karen Lee for the bug report.

  * Bug fix for first-brillouin-zone-k, thanks to Mischa Megens.

  * Bug fix: get-field and compute-field-integral now use fields with the same
    phase as the outputted fields. Thanks to Jim West for the bug report.

  * Miscellaneous bug fixes.

## MPB 1.4.2

3/3/2003.

  * Interactive prompt is now `mpb>` not `guile>`.

  * Output `freqs:` line lists headings as `k1, k2, k3` instead of
    `kx, ky, kz` since they are in reciprocal-lattice, not cartesian,
    coordinates.  Thanks to Theis Peter Hanson for the suggestion.

  * Bug fix in find-k for non-orthogonal lattices; thanks to Suxia
    (Susan) Yang for tracking down this bug.

  * Fixed SunOS problem where k vectors along no-size dimensions failed;
    thanks to Benjamin Cowan for the bug report.

  * Fixed find-k to work for band-min > 1; thanks to M. Povinelli for
    the bug report.

  * Fixed find-k to work for thunk band functions (which take no arguments
    and are called only once instead of per-band).

## MPB 1.4.1

9/16/2002.

  * Fixed NaN in field normalization when basis determinant was negative.
    Thanks to Rumen Iliew for the bug report.

  * Fixed compatibility problems with versions of Guile prior to 1.4.
    Thanks to Cazimir G. Bostan for bug reports.

  * Don't resize lattice basis for grid-size == 1 unless no-size was
    explicitly specified; thanks to Tairan Wang for the suggestion.

## MPB 1.4

9/12/2002.

  * New find-k routine to find k as a function of frequency, instead of
    vice-versa.

  * The Great Field Renormalization: all fields are now normalized
    to have unit *integral* of their energy density (instead of
    unit sum over the grid points), which is much more useful e.g. for
    perturbation theory.  (See field normalization section of manual.)

  * You can now save fields in Scheme variables to perform computations
    combining different fields.  Example routines, e.g. to output the
    Poynting vector, are included.

  * Functions to export (and import) the raw eigenvectors (planewave
    amplitudes), as well as to compute dot products of eigenvectors
    from different k-points (e.g. for detecting band crossings).

  * allow-negative-epsilon function to enable negative-dielectric support.

  * Added examples/dos.scm to compute density of states via simple
    Gaussian histogram, suggested by Xavier Gonze and Doug Allan.

  * Bug fix: allow real offdiagonal epsilon elements without requiring
    --with-hermitian-eps.  Thanks to Doug Allan for the bug report.

  * Eliminated floating-point error on Alpha for homogeneous structure.
    Thanks to F. Lopez-Tejeira for the bug report.

  * Added man page for mpb-split.

  * Minor installation fixes.

## MPB 1.3

3/10/2002.

  * You can now specify the grid size via the resolution input variable,
    instead of via grid-size.  In this case, you make e.g. a 2d simulation
    by creating a lattice with size no-size in one dimension.  The
    old syntax is still supported, but the new style is encouraged
    (all examples have been updated to the new style).

  * New functions to retrieve fields, dielectric functions, etcetera
    at any point, interpolated from the grid if necessary; see the
    get-*-point functions in the reference section.

  * New compute-field-integral function, analogous to compute-energy-integral;
    thanks to Marin Soljacic for the suggestion.

  * Support Scheme complex numbers where appropriate (e.g. in epsilon-offdiag
    or in the new field integration functions).

  * Got rid of NaN when computing the (undefined) group velocity
    for zero-frequency states at the Gamma point; arbitrarily
    return zero here instead.  Thanks to Dmitry N. Chigrin for
    reporting floating-point exceptions on Alphas.

  * Fixed compilation failure for Fortran compilers that use all upper case;
    thanks to Steve Lantz of Cornell.

  * Added "Fun with Fortran" section to installation manual describing
    common Fortran pitfalls; thanks to Steve Lantz for the suggestion.

  * Improved BLAS/LAPACK detection; new --with-blas and --with-lapack
    options to specify these libraries manually.

  * Shortened --with-hermitian-epsilon configure option to 
    --with-hermitian-eps.

  * The data-analysis tutorial is now consistent with h5topng 1.7.

  * Use new API from libctl 2.0.

## MPB 1.2.2

12/7/2001.

  * Fixed bug that caused erroneous/failed convergence when EVEN-Y/ODD-Y
    constraints were used in three dimensions.  Thanks to Rumen Iliew
    for the bug report.

  * Added convenience functions run-yeven, run-yodd, run-yeven-zeven, ...

## MPB 1.2.1

11/20/2001.

  * Fixed serious crashing bug in 1.2; thanks to Karl Koch for the bug report.

## MPB 1.2

11/15/2001.

  * Added new y-parity computation and constraints.  See the new
    run-parity function, which allows you to simultaneously specify
    the parity through the y=0 and z=0 planes, for symmetric structures.
    See also the display-yparities function.

  * z parity is no longer computed by default; see the new
    display-yparities and display-zparities functions to pass to (run).

  * Return more-accurate average epsilon, fill factor, and scalar
    epsilon values (eigenfrequencies are not affected).  Thanks
    to Mischa Megens for bugging me.

  * Now outputs D and H in consistent units (previously, D was multiplied
    by a factor of -frequency).  Thanks to Michelle Povinelli for worrying.

  * epsilon.h5 file now includes extra datasets for all components of the
    effective dielectric tensor.  (This feature is not yet supported if
    you configure --with-inv-symmetry --with-hermitian-epsilon.)

  * run-polarization is replaced by run-parity, and run-even/run-odd are
    deprecated in favor of run-zeven/run-zodd.  run-te/run-tm are now
    equivalent to run-zeven/run-zodd when invoked for 3d systems.

  * Noted new basis-size property of geometry-lattice, from libctl 1.5.
    This makes it easier to use conventional units in the fcc lattice.

  * Group-velocity computation no longer silently invalidates fields
    that have been loaded with get-dfield, etcetera.  Thanks to Marin
    Soljacic for the bug report.

  * The configure script now checks that guile is in the $PATH.  Thanks to
    Bing Li and Giridhar Malalahalli for their bug reports.

  * Rotated the W and K points of the diamond-lattice example so that
    they are oriented similarly to those in the Photonic Crystals book
    by Joannopoulos et al. (eigenfrequencies are not affected).  Thanks
    to Robert Sheldon for pointing out that this was confusing.

  * Added honey-rods.ctl example file: a 2d honeycomb lattice of rods.

  * Added line-defect.ctl example file: a line-defect waveguide in
    a 2d triangular lattice of dielectric rods, formed by a missing
    row of rods.

## MPB 1.1.1

7/4/2001.

  * Fixed bug in H-field output that caused subtly incorrect H-field
    files (only) for 3d problems when NOT using mpbi.

  * Fixed bug that caused mpbi to output incorrect results for 1d
    problems (e.g. outputted dielectric functions with zeros).

  * Changed default eigensolver tolerance from 1e-4 to 1e-7.

  * Added retrieve-gap convenience function to return the gap between
    two specified bands.

  * Fixed typo that prevented compilation of MPI (parallel) version.

  * C compiler flags -O3 are no longer used by default, since they
    don't work with some compilers; most of the performance depends
    upon the BLAS and FFTW anyway.  (Users wishing greater optimization
    can set the CFLAGS environment variable.)  Thanks to
    Giridhar Malalahalli for the bug report.

## MPB 1.1

5/6/2001.

  * Added compute-energy-integral function to make it easier to compute
    arbitrary field-energy integrals for perturbation theory; thanks to
    Marin Soljacic for the suggestion.

  * Fixed bug in output-field routines for the case of a nonzero kz
    component, that caused the fields to be multiplied by an exp(ikx)
    phase with a k in the wrong direction.  Thanks to Jesper Riishede
    for the bug report.

## MPB 1.0

2/23/2001.

  * At long last, support for distributed-memory parallel machines
    with MPI.  The computation time (and memory usage) can often improve
    nearly linearly with the number of processors.  Thanks to
    Clarendon Photonics for funding this work.

  * Also added mpb-split script to parallelize in a simpler way, without
    MPI, on e.g. SMP machines, by dividing up the list of k-points
    among a number of serial mpb processes.

  * Fixed bug in mpbi where artifacts could be introduced in 3d field
    and dielectric-function output files.  (This only affected the
    output files, not the frequency eigenvalues, etcetera.)  Thanks
    to Michelle Povinelli for the bug report.

  * Added new material-function material type, so that you can now
    specify that the dielectric tensor be an arbitrary function of
    position.  Thanks to Peter Bermel for needing this.

  * If MPB is configured with the flag --with-hermitian-epsilon, then
    complex-hermitian dielectric tensors (corresponding to magnetic
    materials, which break time-reversal symmetry) are supported.
    Thanks to Shanhui Fan for pestering me about this.

  * Eliminated output-copies input variable; if you want to visualize
    multiple unit cells, you should use mpb-data.

  * Added new `nothing` material that punches a hole through other
    objects to the background.  (This is distinct from default-material
    when epsilon-input-file is used, or for compute-energy-in-objects.)

  * Fixed inability of MPB 0.13 to run under an old version (1.2) of Guile.

  * Now gives an error if k-point or dielectric tensor is incompatible with
    run-te/run-tm, or if the dielectric tensor is not positive-definite.

  * Default to vendor cc instead of gcc, so that C and Fortran compilers
    are in sync.  (We default to the vendor f77 because it was probably
    used to compile LAPACK/BLAS, and Fortran libraries are picky.)

  * The manual now cites our recent publication on the methods behind MPB.

  * Bug fix in compute-energy-in-object-list for non-orthogonal lattices.

  * Bug fix in combine-band-functions and other functions of band functions,
    which did not handle functions of no arguments ("thunks") correctly
    (crashing with an error message).  Thanks to Michelle Povinelli for
    the bug report.

  * Fixed a floating-point sensitivity bug in mpb-data that could cause
    a crash on the Alpha; thanks to Dominique Caron for the bug report
    and debugging information.

## MPB 0.13

1/7/2001.

  * Can now take advantage of inversion symmetry in the geometry, gaining
    at least a factor of two in speed and a factor of two in memory.
    To use this, you configure MPB with --with-inv-symmetry; the resulting
    executable is installed as `mpbi` and only supports inversion symmetry,
    so you will usually want to install the ordinary MPB as well.

  * Added new eigensolver-block-size input variable, so that MPB can
    optionally solve for only a few bands at a time instead of all
    at once, reducing memory requirements and often increasing speed.

  * Improved handling of the singular (zero-frequency) solutions at
    the Gamma (k=0) point.  This k point should no longer converge
    slowly (or cause additional problems in the targeted eigensolver).

  * Manual updates: please see new referencing suggestions; expanded table
    of contents; we now use more conventional units in diamond/fcc example.

  * You can now pass a "thunk" (function of no arguments) to run, and
    it will be evaluated once per k-point (instead of once per band
    per k-point as for ordinary band functions).

  * compute-field-energy function now also returns the fraction of the
    energy in the various field components.  Thanks to Karl Koch
    for the suggestion.

  * The filename-prefix variable is now read each time an output function
    is called, instead of once per (run), so it can be changed frequently
    if desired.  Thanks to Karl Koch for the suggestion.

  * Added first-brillouin-zone function to transform an arbitrary
    k-point into an equivalent point in the first Brillouin zone.
    Thanks to Payam Rabiei for the suggestion.

  * In mpb-data, the center of the output cell is now always identical
    to the origin of the coordinate system.  Thanks to Michelle
    Povinelli for pointing out this deficiency.	

  * Used improved spherical-quadrature formula in computing the
    effective dielectric tensor in 3d; this should increase accuracy
    somewhat at lower grid resolutions.  Thanks to Doug Allan for
    helpful discussions.

## MPB 0.12

7/9/2000.

  * Added fix-*field-phase functions to allow a deterministic phase
    in the output fields, thanks to a suggestion by Doug Allan.

  * Added group-velocity calculation functions (display-group-velocities,
    etcetera).

  * Added -e x,y,z option to mpb-data so that you can now specify
    an orientation of the output cell (e.g. to make the first axis
    the 111 direction of an fcc crystal).

  * Added (index n) substitute for epsilon property of dielectrics,
    equivalent to (epsilon (* n n)).

  * Documented new libgeom features: cone geometric object, coordinate
    conversion functions (reciprocal->lattice, lattice->cartesian, etc.),
    and vector/matrix rotation.

  * compute-field-energy now returns the total, unnormalized energy
    in the corresponding field; combined with compute-energy-in-objects,
    this makes it easy to do some perturbation theory and related
    calculations.

  * Eigensolver improvements.  Periodic reorthogonalization and
    renormalization to combat some numerical problems. New 
    line-minimization code, included with permission from MINPACK-2
    by Jorge More.

  * Fixed breaking of 90-degree rotational symmetry-breaking by the mesh
    in 2d; thanks to Jim West and Doug Allan of Corning for the bug
    report.  (In general, some symmetry-breaking by the discretization
    seems hard to avoid, however.)

  * Fixed bug in field output routines that could cause crashes
    for grid sizes not a multiple of 4.

  * Bug fix in dielectric function construction for 2d systems: we
    now use the xy plane at z=0 as documented, instead of z=-0.5.

## MPB 0.11

2/12/2000.

  * configure script can now detect and link ATLAS 3.0 accelerated BLAS.

  * Added band-range-data output variable.

  * Running mpb-data multiple times on the same file now replaces
    the results of the previous run, instead of appending -new2,
    -new3, etcetera.

  * Fixed bug in run-even/run-odd that could seriously slow or even
    prevent eigensolver convergence.  Thanks to Payam Rabiei for the
    bug report.

  * Fixed compilation --without-hdf5, or when HDF5 is not found.  Thanks
    to Rajesh Rengarajan for the bug report.

## MPB 0.10

1/28/2000.

  * Added mpb-data utility for post-processing data (e.g. for unskewing
    non-orthogonal lattices).  See the data analysis tutorial or
    man mpb-data for more information.

  * Added new data analysis tutorial to the manual, describing how to
    analyze and visualize the results of two sample calculations.

  * Added support for a new material type, dielectric-anisotropic, so
    that you can specify arbitrary real/symmetric dielectric tensors.

  * Added new output-at-kpoint function to make it easier to output
    fields only at a single k-point in a band-structure calculation.

  * When outputting fields, output all field components (x, y, z, and
    real and imaginary parts) to a single HDF5 file.  Also include info
    on the lattice and k-point vectors to facilitate post-processing.

  * Added new subsection to the installation manual describing some
    generic installation path issues on Unix that were confusing people.

  * Use CPPFLAGS environment variable instead of the less-standard
    INCLUDES to pass -I flags to the configure script (for header files
    in non-standard locations).

  * Added diamond.ctl example file for a 3d diamond (fcc) lattice of spheres.

  * Added (brief) mpb man page.

  * Fixed z-parity output and run-even/odd functions for 2d grids.

  * Fixed bug in output-dpwr-in-objects.  Thanks to Mihai Ibanescu for
    the bug report.

  * Compilation fixes.  We need to set SHELL in the Makefile for make on
    some systems.  Also added rule to insure ctl-io.h is created before
    main.c is compiled.  Thanks to Christoph Becher for the bug reports.

## MPB 0.9.1

1/7/2000.

  * Fixed eigensolver bug where special handling of Gamma (k=0) point could
    screw up convergence for subsequent k-points, causing incorrect results.

  * Fixed behavior of filename-prefix input variable; thanks to Karl Koch
    for the bug report.

## MPB 0.9

1/2/2000.

  * Added run-even and run-odd functions, so you can now compute only
    even/odd states (with respect to a z=0 mirror plane) in systems with
    sufficient symmetry.  See also the new z-parity output variable.

  * Added epsilon-input-file variable, so that you can now read an
    arbitrary dielectric function from a file.

  * Field file names now include the polarization (e.g. `.tm`).

  * Some optimizations in the eigensolver.

  * Some documentation improvements; thanks to Edmond Chow for his comments.

  * configure should work even when there is no Fortran compiler on your
    system (assuming your BLAS, etc., libraries work without Fortran libs).
    Thanks to Antti Renko for the bug report.

  * Fixed problems detecting BLAS and LAPACK shared libraries in configure.
    Thanks to Karri Varris for the bug report.

  * Fixed trailing spaces in sed command, which were breaking `make install`
    on some systems.  Thanks to Ron Chase for the bug report.

## MPB 0.8.1

11/22/1999.

  * Added output-hfield-x, output-dfield-y, etcetera, functions for
    outputting only specific field components (see manual reference section).

  * Sped up HDF5 field output routines.

  * Added output-copies variable to set the number of periods output by
    the band output functions (see manual reference section).

## MPB 0.8

11/19/1999.

  * Initial public release.
