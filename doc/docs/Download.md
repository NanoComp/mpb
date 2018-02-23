---
# Download
---

The latest development source is on [GitHub](https://github.com/stevengj/mpb).

The current stable release is **version 1.6.1** which can be downloaded from:

- <https://github.com/stevengj/mpb/releases/download/v1.6.1/mpb-1.6.1.tar.gz>

Refer to [NEWS](https://github.com/stevengj/mpb/blob/master/NEWS.md) to see what's new in this version, and be sure to read the [Installation](Installation.md) for how to compile and install it.

Please subscribe to the **mpb-announce** mailing list to receive notifications when new versions are released:

-   [mpb-announce mailing list](http://ab-initio.mit.edu/cgi-bin/mailman/listinfo/mpb-announce)

MPB on Amazon Web Services (AWS)
---------------------------------

The most recent, stable version of MPB preinstalled on [Ubuntu](https://en.wikipedia.org/wiki/Ubuntu) can be accessed for free on Amazon Web Services (AWS) Elastic Compute Cloud (EC2) as an [Amazon Machine Image (AMI)](https://aws.amazon.com/marketplace/pp/B01KHWH0AS) provided by [Simpetus](http://www.simpetus.com/launchsims.html).

Precompiled MPB Packages for and Ubuntu
---------------------------------------

Precompiled packages of MPB are available for [Ubuntu](https://en.wikipedia.org/wiki/Ubuntu) as **mpb**. The Ubuntu 16.10 package is available in the [science](https://packages.ubuntu.com/yakkety/mpb) repository. We highly recommend using Ubuntu as MPB and all of its dependencies can be installed using just one line:

```
sudo apt-get install mpb h5utils
```

You can also install the [parallel version of MPB](https://packages.ubuntu.com/trusty/science/mpb-mpi) which is based on [MPICH](https://www.mpich.org/) using:

```
sudo apt-get install mpb-mpi
```