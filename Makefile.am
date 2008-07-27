ACLOCAL_AMFLAGS=-I ./m4

SUBDIRS = src tests mpb utils # meb

EXTRA_DIST = COPYRIGHT TODO doc examples m4

dist-hook:
	darcs changes --summary > $(distdir)/ChangeLog
