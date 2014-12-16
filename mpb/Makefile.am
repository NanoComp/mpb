# libctl: flexible Guile-based control files for scientific software
# Copyright (C) 1998, 1999-2012, Steven G. Johnson
#
# This file may be used without restriction.  It is in the public
# domain, and is NOT restricted by the terms of any GNU license.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# Steven G. Johnson can be contacted at stevenj@alum.mit.edu.

##############################################################################

bin_PROGRAMS = mpb@MPB_SUFFIX@

SPECIFICATION_FILE = mpb.scm
EXTRA_DIST = mpb.scm.in epsilon.c mu.c

nodist_pkgdata_DATA = $(SPECIFICATION_FILE)

MY_SOURCES = medium.c epsilon_file.c field-smob.c fields.c	\
material_grid.c material_grid_opt.c matrix-smob.c mpb.c field-smob.h matrix-smob.h mpb.h my-smob.h

MY_LIBS = $(top_builddir)/src/matrixio/libmatrixio.a $(top_builddir)/src/libmpb@MPB_SUFFIX@.la $(NLOPT_LIB)
MY_CPPFLAGS = -I$(top_srcdir)/src/util -I$(top_srcdir)/src/matrices -I$(top_srcdir)/src/matrixio -I$(top_srcdir)/src/maxwell
MY_DEFS = -DHAVE_CTL_HOOKS=1 -DHAVE_CTL_EXPORT_HOOK=1

# what is printed out when invoking your program with --version:
VERSION_STRING = "mpb@MPB_SUFFIX@ @VERSION@, Copyright (C) 1999-2012, MIT"

check-local: mpb@MPB_SUFFIX@
	./mpb@MPB_SUFFIX@ $(top_srcdir)/examples/check.ctl
	./mpb@MPB_SUFFIX@ force-mu?=true $(top_srcdir)/examples/check.ctl

dist_man_MANS = mpb-split.1 mpb.1

if !MPI
bin_SCRIPTS = mpb@MPB_SUFFIX@-split
endif

mpb@MPB_SUFFIX@-split: mpb-split-preinstall
	sed "s%\./mpb%"$(prefix)/bin/"mpb@MPB_SUFFIX@%g" \
            mpb-split-preinstall > $@
	chmod u+x $@

# The following variables should be detected and set by autoconf:

# libctl install. dir., e.g. /usr/local/share/libctl
LIBCTL_DIR = @LIBCTL_DIR@

# gen-ctl-io program
GEN_CTL_IO = @GEN_CTL_IO@

##############################################################################
#   don't (normally) edit below except to replace "mpb" with your name   #
##############################################################################

CTL_DEFS = -DCTL_SCM='"'$(LIBCTL_DIR)/base/ctl.scm'"' \
       -DINCLUDE_SCM='"'$(LIBCTL_DIR)/base/include.scm'"' \
          -DSPEC_SCM='"'$(pkgdatadir)/$(SPECIFICATION_FILE)'"' \
    -DVERSION_STRING='"'$(VERSION_STRING)'"'

mpb@MPB_SUFFIX@_SOURCES = $(MY_SOURCES)
nodist_mpb@MPB_SUFFIX@_SOURCES = main.c geom.c ctl-io.h ctl-io.c
BUILT_SOURCES = $(nodist_mpb@MPB_SUFFIX@_SOURCES) mu.c
mpb@MPB_SUFFIX@_LDADD = $(MY_LIBS)
mpb@MPB_SUFFIX@_LDFLAGS = $(MY_LDFLAGS)
mpb@MPB_SUFFIX@_CPPFLAGS = $(MY_CPPFLAGS) $(MY_DEFS) $(CTL_DEFS)  $(CTL_H_CPPFLAG)

if MAINTAINER_MODE

mu.c: epsilon.c
	$(PERL) -pe 's/epsilon/mu/g;s/eps_inv/mu_inv/g' $< > $@

endif

main.c: $(LIBCTL_DIR)/base/main.c
	cp -f $(LIBCTL_DIR)/base/main.c $@

geom.c: $(LIBCTL_DIR)/utils/geom.c
	cp -f $(LIBCTL_DIR)/utils/geom.c $@

ctl-io.c: $(SPECIFICATION_FILE) $(LIBCTL_DIR)/utils/geom.scm
	$(GEN_CTL_IO) --code -o $@ $(SPECIFICATION_FILE) $(LIBCTL_DIR)

ctl-io.h: $(SPECIFICATION_FILE) $(LIBCTL_DIR)/utils/geom.scm
	$(GEN_CTL_IO) --header -o $@ $(SPECIFICATION_FILE) $(LIBCTL_DIR)

clean-local:
	rm -f ctl-io.* main.* geom.*
