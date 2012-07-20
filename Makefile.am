ACLOCAL_AMFLAGS=-I m4

SUBDIRS = src tests mpb utils examples # meb

EXTRA_DIST = COPYRIGHT TODO doc m4

dist-hook:
	darcs changes --summary > $(distdir)/ChangeLog
