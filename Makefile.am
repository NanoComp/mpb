ACLOCAL_AMFLAGS=-I m4

SUBDIRS = src tests # meb

if WITH_LIBCTL
SUBDIRS += mpb examples
endif
if WITH_LIBCTLGEOM
SUBDIRS += utils
endif

EXTRA_DIST = COPYRIGHT TODO m4 README.md

dist-hook:
	git log > $(distdir)/ChangeLog
