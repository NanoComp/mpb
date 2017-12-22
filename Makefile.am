ACLOCAL_AMFLAGS=-I m4

SUBDIRS = src tests utils # meb

if WITH_LIBCTL
SUBDIRS += mpb examples
endif

EXTRA_DIST = COPYRIGHT TODO m4 README.md

dist-hook:
	git log > $(distdir)/ChangeLog
