---
# MPB
---

The **MIT Photonic-Bands** (**MPB**) package is a [free](http://www.gnu.org/philosophy/free-sw.en.html) program for computing the band structures, or dispersion relations, and electromagnetic modes of periodic dielectric structures, on both serial and parallel computers. MPB computes definite-frequency eigenstates, or harmonic modes, of [Maxwell's equations](https://en.wikipedia.org/wiki/Maxwell%27s_equations) in periodic dielectric structures for arbitrary wavevectors, using fully-vectorial and three-dimensional methods. It is applicable to many problems in optics, such as waveguides and resonator systems, and [photonic crystals](http://ab-initio.mit.edu/book). For example, it can solve for the modes of waveguides with arbitrary cross-sections.

See also our complementary [Meep](http://meep.readthedocs.io/en/latest/Meep/) package for time-domain simulations, reflection/transmission spectra, etc.

Features
--------

-   Free software under the GNU General Public License. See the [License and Copyright](License_and_Copyright.md).
-   Fully-vectorial, three-dimensional calculation. Iterative eigensolver techniques are employed to make large, three-dimensional calculations possible. Can also handle 2D and 1D problems.
-   Direct, frequency-domain eigensolver as opposed to indirect methods, e.g. time-domain. For one thing, this means that you get both eigenvalues (frequencies) and eigenstates (electromagnetic modes) at the same time. See a [comparison of time-domain and frequency-domain techniques](Introduction.md#frequency-domain-vs-time-domain).
-   Targeted eigensolver. Normally, iterative eigensolvers provide you with the states (optical bands/modes) with the lowest few frequencies. Our software can alternatively compute the modes whose frequencies are closest to a specified target frequency. This greatly reduces the number of bands that must be computed in guided or resonant mode calculations.
-   Flexible, scriptable user interface based on the [Guile](http://www.gnu.org/software/guile/) extension & scripting language.
-   Support for arbitrary, anisotropic dielectric structures including gyrotropic/magneto-optic materials and non-orthogonal unit cells.
-   Field output in [HDF5](https://support.hdfgroup.org/HDF5/) format for input into many popular graphing and visualization tools.
-   Portable to most Unix-like operating systems. See the [installation guide](Installation.md).
-   Support for parallel machines with MPI.

To give you some feel for how long these calculations take, let us consider one typical data point. For the 3d band-structure of a [diamond lattice of dielectric spheres in air](Data_Analysis_Tutorial.md#diamond-lattice-of-spheres), computing the lowest 10 bands on a 16×16×16 grid at 31 k-points, MPB took 8 seconds on a 2.8 GHz AMD Opteron under Linux with the [ATLAS](http://www.netlib.org/atlas/) optimized BLAS library. Thus, at each k-point, MPB was minimizing a function with 81920 degrees of freedom in 0.26 seconds on average.

MPB Download
------------

You can [download](Download.md) the full source code in ANSI C for MPB under the [GPL](License_and_Copyright.md). The [installation guide](Installation.md) describes how to install it; mainly, this consists of downloading and installing various prerequisites if you do not have them already.

The current version is **1.5**. See the [release notes](Release_Notes.md) for what is new in each version.

You can also download the latest development sources from [MPB on Github](https://github.com/stevengj/mpb).

Documentation
-------------

This manual is readable online and is also part of the code repository. The [tutorials](Scheme_Tutorial.md) demonstrate what it is like to use the program. You may be also interested in the [libctl manual](http://ab-initio.mit.edu/libctl), which describes a Guile/Scheme-based scripting library that we build our interface on top of, and also additional information on [Guile and Scheme](Guile_and_Scheme_Information.md).

We have published a paper on the computational methods underlying MPB:

Steven G. Johnson and J. D. Joannopoulos, [Block-iterative frequency-domain methods for Maxwell's equations in a planewave basis](http://www.opticsinfobase.org/abstract.cfm?URI=oe-8-3-173), *Optics Express* **8**, no. 3, 173-190 (2001).

### Mailing Lists

The MPB mailing lists and their archives are another source of information about MPB.

Subscribe to the read-only [mpb-announce mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/mpb-announce) to receive notifications of updates and releases. Subscribe to the unmoderated [mpb-discuss mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/mpb-discuss) for discussions about using MPB. Archives are available [here](http://www.mail-archive.com/mpb-discuss@ab-initio.mit.edu/). You can also read and post to the list via the [gmane.comp.science.photonic-bands](news://news.gmane.org/gmane.comp.science.photonic-bands) newsgroup from [Gmane](http://www.gmane.org/).

### Bug Reports and Feature Requests

For bug reports and feature requests, please [file an MPB Github issue](https://github.com/stevengj/mpb/issues).

Acknowledgements
----------------

Many people and groups have contributed to the development of this software, both directly and indirectly. Please see the [acknowledgements section](Acknowledgements.md) of the manual for those to whom we feel especially grateful.

Contacts and Feedback
---------------------

If you have questions or problems regarding MPB, you are encouraged to query the [mailing list](https://www.mail-archive.com/mpb-discuss@ab-initio.mit.edu/).

For professional consulting as well as free access to MPB in the public cloud via Amazon Web Services (AWS), see [Simpetus](http://www.simpetuscloud.com).