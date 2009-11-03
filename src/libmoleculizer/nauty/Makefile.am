## Process this file with automake to produce Makefile.in.
noinst_LTLIBRARIES = libmoleculizer_nauty.la

VERSION = @INSTALLATION_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_nautydir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/nauty

libmoleculizer_nauty_la_SOURCES =\
nauty.c \
nautil.c \
naugraph.c

libmoleculizer_nauty_HEADERS = \
gtools.h \
gutils.h \
naugroup.h \
nausparse.h \
nautaux.h \
nautinv.h \
naututil.h \
nauty.h \
planarity.h \
rng.h