[![Latest Docs](https://readthedocs.org/projects/pip/badge/?version=latest)](http://mpb.readthedocs.io/en/latest/)
[![Build Status](https://travis-ci.org/stevengj/mpb.svg?branch=master)](https://travis-ci.org/stevengj/mpb)

MPB is a free and open-source software package for computing electromagnetic band structures and modes.

**Features**

-   **Free and open-source software** under the [GNU GPL](https://en.wikipedia.org/wiki/GNU_General_Public_License).
-   Complete **scriptability** via [Python](Python_Tutorial) or [Scheme](Scheme_User_Interface) APIs.
-   Portable to any Unix-like system such as [Linux](https://en.wikipedia.org/wiki/Linux), [macOS](https://en.wikipedia.org/wiki/MacOS), and [FreeBSD](https://en.wikipedia.org/wiki/FreeBSD).
-   Distributed memory **parallelism** on any system supporting the [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface) standard.
-   Fully-vectorial **1d, 2d, 3d** calculations. Iterative eigensolver techniques are employed to make large calculations possible.
-   **Direct, frequency-domain eigensolver** as opposed to indirect methods, e.g. time-domain. This means that you get both eigenvalues (frequencies) and eigenstates (electromagnetic modes) at the same time. See [comparison of time-domain and frequency-domain techniques](Introduction.md#frequency-domain-vs-time-domain).
-   **Targeted eigensolver**. Iterative eigensolvers normally compute states (harmonic modes) with the lowest few frequencies. MPB can alternatively compute the modes whose frequencies are closest to a specified target frequency. This greatly reduces the number of bands that must be computed in guided or resonant mode calculations.
-   Support for arbitrary, **anisotropic** dielectrics including **gyrotropic/magneto-optic** materials and **non-orthogonal** unit cells.
-   Field output in the [HDF5](https://support.hdfgroup.org/HDF5/) data format.

# Documentation

See the [manual on readthedocs](https://mpb.readthedocs.io/en/latest) for the latest documentation.