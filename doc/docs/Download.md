---
# Download
---

The latest development sources are available on [GitHub](https://github.com/stevengj/mpb).

The latest release of MPB is **version 1.5** which can be downloaded from:

- <http://ab-initio.mit.edu/mpb/mpb-1.5.tar.gz>

Older releases may be found at <http://ab-initio.mit.edu/mpb/old>

MPB is free software under the [GNU GPL](License_and_Copyright.md).

Refer to the [Release Notes](Release_Notes.md) to see what's new in this version, and be sure to read the [installation guide](Installation.md) for how to compile and install it.

Please subscribe to the **mpb-announce** mailing list to receive notifications when new versions are released:

-   [mpb-announce mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/mpb-announce)

MPB on Amazon Web Services (AWS)
---------------------------------

The latest version of Meep preinstalled on [Ubuntu](https://en.wikipedia.org/wiki/Ubuntu) can be accessed for free on Amazon Web Services (AWS) Elastic Compute Cloud (EC2) as an [Amazon Machine Image (AMI)](https://aws.amazon.com/marketplace/pp/B01KHWH0AS). Instructions for accessing this AMI are provided by [Simpetus](http://www.simpetuscloud.com/launchsims.html).

Precompiled MPB Packages for Debian and Ubuntu
-----------------------------------------------

Precompiled packages of MPB are available for [Debian](https://en.wikipedia.org/wiki/Debian) and [Ubuntu](https://en.wikipedia.org/wiki/Ubuntu) as **mpb**. The Debian package is available in the [testing](http://packages.debian.org/testing/science/mpb) repository. We highly recommend using Debian or Ubuntu, as in Debian or Ubuntu the Meep software and all of its dependencies can be installed using just one line:

```
sudo apt-get install mpb h5utils
```

You can also install the [parallel version of MPB](http://packages.debian.org/testing/science/mpb-mpi) which is based on [MPICH](https://www.mpich.org/) using:

```
sudo apt-get install mpb-mpi
```

Python User Interface
----------------

A Python user interface for Meep is currently under development by [Simpetus](http://www.simpetuscloud.com).