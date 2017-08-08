---
title: MPB on Debian
permalink: /MPB_on_Debian/
---

The most tedious and error-prone part of MPB is installing all of the [prerequisite packages](/MPB_Installation "wikilink"). For users of [Debian GNU/Linux](http://www.debian.org/), however, this process can be entirely automated.

[Josselin Mouette](mailto:josselin.mouette@ens-lyon.org) has graciously volunteered to provide Debian packages of MPB. These packages are now in the official Debian:

-   [mpb](http://packages.debian.org/testing/science/mpb.html)
-   [mpb-doc](http://packages.debian.org/testing/doc/mpb-doc.html)

You can therefore simply do:

`apt-get install mpb`

as `root` to install MPB and all its prerequisites.

You can also

`apt-get install h5utils`

to install the [h5utils](/h5utils "wikilink") package, which is useful in conjunction with MPB.

[Category:MPB](/Category:MPB "wikilink")