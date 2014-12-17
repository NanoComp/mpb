ACLOCAL_AMFLAGS=-I m4

SUBDIRS = src tests mpb utils examples # meb

EXTRA_DIST = COPYRIGHT TODO m4 README.md

dist-hook:
	git log > $(distdir)/ChangeLog
