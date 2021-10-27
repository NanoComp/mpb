---
# MPB
---

**MPB** is a free and open-source software package for computing the band structures, or dispersion relations, and electromagnetic modes of periodic dielectric structures, on both serial and parallel computers. MPB is an acronym for *MIT Photonic Bands*. MPB computes definite-frequency eigenstates, or harmonic modes, of [Maxwell's equations](https://en.wikipedia.org/wiki/Maxwell%27s_equations) for periodic dielectric structures of lossless, wavelength-independent, anisotropic $\varepsilon$ and $\mu$ for arbitrary wavevectors, using fully-vectorial and three-dimensional methods. It is applicable to many problems in optics, such as waveguides and resonator systems, and [photonic crystals](http://ab-initio.mit.edu/book).

See also the complementary [Meep](https://meep.readthedocs.io/) package for time-domain simulations, reflection/transmission spectra, etc.

Features
--------

-   **Free and open-source software** under the [GNU GPL](https://en.wikipedia.org/wiki/GNU_General_Public_License).
-   Complete **scriptability** via [Python](Python_Tutorial) or [Scheme](Scheme_User_Interface) APIs.
-   Portable to any Unix-like system such as [Linux](https://en.wikipedia.org/wiki/Linux), [macOS](https://en.wikipedia.org/wiki/MacOS), and [FreeBSD](https://en.wikipedia.org/wiki/FreeBSD).
-   Distributed memory **parallelism** on any system supporting the [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface) standard.
-   Fully-vectorial **1d, 2d, 3d** calculations. Iterative eigensolver techniques are employed to make large calculations possible.
-   **Direct, frequency-domain eigensolver** as opposed to indirect methods, e.g. time-domain. This means that you get both eigenvalues (frequencies) and eigenstates (electromagnetic modes) at the same time. See [comparison of time-domain and frequency-domain techniques](Introduction.md#frequency-domain-vs-time-domain).
-   **Targeted eigensolver**. Iterative eigensolvers normally compute states (harmonic modes) with the lowest few frequencies. MPB can alternatively compute the modes whose frequencies are closest to a specified target frequency. This greatly reduces the number of bands that must be computed in guided or resonant mode calculations.
-   Support for arbitrary, **anisotropic** dielectrics including **gyrotropic/magneto-optic** materials and **non-orthogonal** unit cells. Lossy and wavelength-dependent $\varepsilon$ and $\mu$ are not supported.
-   Field output in the [HDF5](https://support.hdfgroup.org/HDF5/) data format.

To give you some feel for how long these calculations take, let us consider one typical data point. For the 3d band-structure of a [diamond lattice of dielectric spheres in air](Data_Analysis_Tutorial.md#diamond-lattice-of-spheres), computing the lowest 10 bands on a 16×16×16 grid at 31 k-points, MPB took 8 seconds on a 2.8 GHz AMD Opteron under Debian with the [ATLAS](http://www.netlib.org/atlas/) optimized BLAS library. Thus, at each k-point, MPB was minimizing a function with 81920 degrees of freedom in 0.26 seconds on average.

Download
------------

The development repository is on on [GitHub](https://github.com/NanoComp/mpb). Gzipped tarballs of stable versions are available in [Download](Download.md). The release history is described in [NEWS](https://github.com/NanoComp/mpb/blob/master/NEWS.md). Installation instructions are in [Installation](Installation.md).

Documentation
-------------

See the navigation sidebar at left. In particular, the [Introduction](Introduction) and [Tutorial](Scheme_Tutorial.md) are the most important things to review.

Please [cite the reference publication](Acknowledgements.md#referencing) in any publication for which you found MPB useful.

### Mailing Lists

Subscribe to the read-only [mpb-announce mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/mpb-announce) to receive notifications of updates and releases. Subscribe to the [mpb-discuss mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/mpb-discuss) for discussions regarding MPB. The [mpb-discuss archives](http://www.mail-archive.com/mpb-discuss@ab-initio.mit.edu/) includes all postings since 2006 spanning a large number and variety of discussion topics related to installation, setting up simulations, post-processing output, etc.

### Bug Reports and Feature Requests

For bug reports and feature requests, please file a [GitHub issue](https://github.com/NanoComp/mpb/issues).

Acknowledgements
----------------

The MPB project is maintained by [Simpetus](http://www.simpetus.com) and the developer community on [GitHub](https://github.com/NanoComp/mpb). [Acknowledgements](Acknowledgements.md) provides a complete listing of the project contributors.

Contacts and Feedback
---------------------

If you have questions or problems regarding MPB, you are encouraged to query the [mailing list](https://www.mail-archive.com/mpb-discuss@ab-initio.mit.edu/).

Professional consulting services for photonic design and modeling including development of custom, turn-key simulation modules, training, technical support, and access to MPB in the public cloud via Amazon Web Services (AWS) are provided by [Simpetus](http://www.simpetus.com).
