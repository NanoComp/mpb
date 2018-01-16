ACLOCAL_AMFLAGS=-I m4

SUBDIRS = src tests mpb utils examples

EXTRA_DIST = COPYRIGHT TODO m4 README.md NEWS.md

dist-hook:
	git log > $(distdir)/ChangeLog
