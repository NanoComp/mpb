ACLOCAL_AMFLAGS=-I m4

SUBDIRS = src tests

if WITH_LIBCTL
SUBDIRS += mpb examples
endif
if WITH_LIBCTLGEOM
SUBDIRS += utils
endif

EXTRA_DIST = COPYRIGHT TODO m4 README.md NEWS.md

# make dist complains about missing NEWS file
NEWS: NEWS.md
	cp NEWS.md $@
