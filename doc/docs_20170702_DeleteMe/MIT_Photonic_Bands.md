---
title: MIT Photonic Bands
permalink: /MIT_Photonic_Bands/
---

[500px|center](/Image:Mpb-logo.jpg "wikilink")

The **MIT Photonic-Bands** (**MPB**) package is a [free](http://www.gnu.org/philosophy/free-sw.html) program for computing the band structures (dispersion relations) and electromagnetic modes of periodic dielectric structures, on both serial and parallel computers. It was developed by [Steven G. Johnson](http://math.mit.edu/~stevenj) at [MIT](http://web.mit.edu/) along with the Joannopoulos [Ab Initio Physics](http://ab-initio.mit.edu/) group.

This program computes definite-frequency eigenstates (harmonic modes) of [Maxwell's equations](/w:Maxwell's_equations "wikilink") in periodic dielectric structures for arbitrary wavevectors, using fully-vectorial and three-dimensional methods. It is especially designed for the study of [photonic crystals](http://ab-initio.mit.edu/book/) (a.k.a. photonic band-gap materials), but is also applicable to many other problems in optics, such as waveguides and resonator systems. (For example, it can solve for the modes of waveguides with arbitrary cross-sections.)

See also our complementary [Meep](/Meep "wikilink") package for time-domain simulations, reflection/transmission spectra, etc.

[MPB download](/MPB_download "wikilink")
----------------------------------------

You can [download](/MPB_download "wikilink") the full source code (in ANSI C) for MPB under the [free GPL license](/MPB_License_and_Copyright "wikilink"). The [installation](/MPB_Installation "wikilink") section of the manual describes how to install it; mainly, this consists of downloading and installing various prerequisites if you do not have them already.

### Latest version: 1.5

The [MPB release notes](/MPB_release_notes "wikilink") describe what is new in each version; the current version of MPB is 1.5.

You can also download the latest development sources from [MPB on Github](https://github.com/stevengj/mpb).

Documentation
-------------

The [MPB manual](/MPB_manual "wikilink") is readable online. Note especially the [tutorial](/MPB_User_Tutorial "wikilink") section of the manual, to get a flavor of what it is like to use the program. You may be also interested in the [libctl manual](/libctl_manual "wikilink"), which describes a Guile/Scheme-based scripting library that we build our interface on top of, and also a collection of [Guile and Scheme links](/Guile_and_Scheme_links "wikilink").

We have published a paper (available online) on the computational methods underlying MPB:


Steven G. Johnson and J. D. Joannopoulos, "[Block-iterative frequency-domain methods for Maxwell's equations in a planewave basis](http://www.opticsinfobase.org/abstract.cfm?URI=oe-8-3-173)," *Optics Express* **8**, no. 3, 173-190 (2001).

See also our [referencing suggestions](/Citing_MPB "wikilink") for how to cite MPB in your work.

### Mailing Lists

The MPB mailing lists and their archives are another source of information about MPB.

Subscribe to the read-only [mpb-announce mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/mpb-announce) to receive notifications of updates and releases. Subscribe to the unmoderated [mpb-discuss mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/mpb-discuss) for discussions about using MPB. Archives are available [here](http://www.mail-archive.com/mpb-discuss@ab-initio.mit.edu/). You can also read and post to the list via the [gmane.comp.science.photonic-bands](news://news.gmane.org/gmane.comp.science.photonic-bands) newsgroup from [Gmane](http://www.gmane.org/).

### Bug reports and feature requests

For bug reports and feature requests, please [file an MPB Github issue](https://github.com/stevengj/mpb/issues).

Features
--------

Some of the more noteworthy features of the MIT Photonic-Bands package are:

-   Fully-vectorial, three-dimensional calculation. Iterative eigensolver techniques are employed to make large, three-dimensional calculations possible. (Can handle 2D and 1D problems too, of course.)
-   Direct, frequency-domain eigensolver (as opposed to indirect methods, e.g. time-domain). For one thing, this means that you get both eigenvalues (frequencies) and eigenstates (electromagnetic modes) at the same time. (See also a [comparison of time-domain and frequency-domain techniques](/MPB_Introduction#Frequency-Domain_vs._Time-Domain "wikilink") in the [MPB manual](/MPB_manual "wikilink").)
-   Targeted eigensolver. Normally, iterative eigensolvers provide you with the states (optical bands/modes) with the lowest few frequencies. Our software can alternatively compute the modes whose frequencies are closest to a specified target frequency. This greatly reduces the number of bands that must be computed in guided or resonant mode calculations.
-   Flexible, scriptable user interface based upon the [GNU Guile](http://www.gnu.org/software/guile/) extension & scripting language.
-   Support for arbitrary, anisotropic dielectric structures (including gyrotropic/magneto-optic materials) and non-orthogonal unit cells.
-   Field output in [HDF](http://hdf.ncsa.uiuc.edu/) format for input into many popular graphing and visualization tools.
-   Portable to most any Unix-like operating system; tested under Linux, AIX, IRIX, and Tru64 (*née* Digital) Unix. See also the [installation section](/MPB_Installation "wikilink") of the manual.
-   Support for parallel machines with MPI. (Tested on an SGI [Origin2000](http://scv.bu.edu/SCV/Origin2000/) and on an SMP Linux machine with [MPICH](http://www-unix.mcs.anl.gov/mpi/mpich/).)
-   Free software under the GNU General Public License. (See also the [MPB License and Copyright](/MPB_License_and_Copyright "wikilink") section of the manual.)

To give you some feel for how long these calculations take, let us consider one typical data point. For the 3d band-structure of a [diamond lattice of dielectric spheres in air](/MPB_Data_Analysis_Tutorial#Diamond_Lattice_of_Spheres "wikilink"), computing the lowest 10 bands on a 16×16×16 grid at 31 k-points, MPB 1.0 took 2 minutes on a 550MHz Pentium-III under Linux with the [ATLAS](http://www.netlib.org/atlas/) optimized BLAS library. (Thus, at each k-point, MPB was minimizing a function with 81920 degrees of freedom in 4 seconds on average.)

[Acknowledgements](/MPB_Acknowledgements "wikilink")
----------------------------------------------------

Many people and groups have contributed to the development of this software, both directly and indirectly. Please see the [acknowledgements section](/MPB_Acknowledgements "wikilink") of the manual for those to whom we feel especially grateful.

Contacts and Feedback
---------------------

If you have questions or problems regarding MPB, you are encouraged to query the [mailing list](https://www.mail-archive.com/mpb-discuss@ab-initio.mit.edu/).

For professional consulting as well as free access to MPB in the public cloud via Amazon Web Services (AWS), see [Simpetus](http://www.simpetuscloud.com).

[Category:MPB](/Category:MPB "wikilink")