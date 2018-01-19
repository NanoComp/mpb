---
# MPB
---

**MPB** is a free software package for computing the band structures, or dispersion relations, and electromagnetic modes of periodic dielectric structures, on both serial and parallel computers. MPB is an acronym for *MIT Photonic Bands*. MPB computes definite-frequency eigenstates, or harmonic modes, of [Maxwell's equations](https://en.wikipedia.org/wiki/Maxwell%27s_equations) in periodic dielectric structures for arbitrary wavevectors, using fully-vectorial and three-dimensional methods. It is applicable to many problems in optics, such as waveguides and resonator systems, and [photonic crystals](http://ab-initio.mit.edu/book). For example, it can solve for the modes of waveguides with arbitrary cross-sections.

See also our complementary [Meep](https://meep.readthedocs.io/) package for time-domain simulations, reflection/transmission spectra, etc.

Features
--------

-   **Free software** under the [GNU GPL](https://en.wikipedia.org/wiki/GNU_General_Public_License).
-   Portable to any Unix-like system such as [Linux](https://en.wikipedia.org/wiki/Linux) and [macOS](https://en.wikipedia.org/wiki/MacOS).
-   Support for **parallel** machines with MPI.
-   Fully-vectorial **1d, 2d, 3d** calculations. Iterative eigensolver techniques are employed to make large calculations possible.
-   **Direct, frequency-domain eigensolver** as opposed to indirect methods, e.g. time-domain. This means that you get both eigenvalues (frequencies) and eigenstates (electromagnetic modes) at the same time. See a [comparison of time-domain and frequency-domain techniques](Introduction.md#frequency-domain-vs-time-domain).
-   **Targeted eigensolver**. Normally, iterative eigensolvers provide you with the states (photonic bands/modes) with the lowest few frequencies. MPB can alternatively compute the modes whose frequencies are closest to a specified target frequency. This greatly reduces the number of bands that must be computed in guided or resonant mode calculations.
-   Flexible, scriptable user interface based on [Scheme](https://en.wikipedia.org/wiki/Scheme_programming_language). A [Python](https://en.wikipedia.org/wiki/Python_programming_language) interface is under development.
-   Support for arbitrary, **anisotropic** dielectricstructures including **gyrotropic/magneto-optic** materials and **non-orthogonal** unit cells.
-   Field output in [HDF5](https://support.hdfgroup.org/HDF5/) format supported by many visualization tools.

To give you some feel for how long these calculations take, let us consider one typical data point. For the 3d band-structure of a [diamond lattice of dielectric spheres in air](Data_Analysis_Tutorial.md#diamond-lattice-of-spheres), computing the lowest 10 bands on a 16×16×16 grid at 31 k-points, MPB took 8 seconds on a 2.8 GHz AMD Opteron under Debian with the [ATLAS](http://www.netlib.org/atlas/) optimized BLAS library. Thus, at each k-point, MPB was minimizing a function with 81920 degrees of freedom in 0.26 seconds on average.

Download
------------

The latest development sources are available on [GitHub](https://github.com/stevengj/mpb). The source tarballs are available on the [Download](Download.md) page. The release history is described in the [NEWS file](https://github.com/stevengj/mpb/blob/master/NEWS.md). The installation instructions can be found in the [Installation](Installation.md) page.

Documentation
-------------

See the navigation sidebar at left. The [Tutorial](Scheme_Tutorial.md) demonstrates what it is like to use the program.

Please [cite MPB](Acknowledgements.md#referencing) in any publication for which you found it useful.

### Mailing Lists

The MPB mailing lists and their archives are another source of information about MPB.

Subscribe to the read-only [mpb-announce mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/mpb-announce) to receive notifications of updates and releases. Subscribe to the unmoderated [mpb-discuss mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/mpb-discuss) for discussions about using MPB. Archives are available [here](http://www.mail-archive.com/mpb-discuss@ab-initio.mit.edu/). You can also read and post to the list via the [gmane.comp.science.photonic-bands](news://news.gmane.org/gmane.comp.science.photonic-bands) newsgroup from [Gmane](http://www.gmane.org/).

### Bug Reports and Feature Requests

For bug reports and feature requests, please [file an MPB Github issue](https://github.com/stevengj/mpb/issues).

Acknowledgements
----------------

The MPB project is maintained by [Simpetus](http://www.simpetuscloud.com) and the open-source community on [GitHub](https://github.com/stevengj/mpb). Please see the [Acknowledgements](Acknowledgements.md) for a more complete listing of those to whom we are grateful.

Contacts and Feedback
---------------------

If you have questions or problems regarding MPB, you are encouraged to query the [mailing list](https://www.mail-archive.com/mpb-discuss@ab-initio.mit.edu/).

Professional consulting services as well as free access to MPB in the public cloud via Amazon Web Services (AWS) are provided by [Simpetus](http://www.simpetuscloud.com).
