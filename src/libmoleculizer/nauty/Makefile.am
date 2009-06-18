## Process this file with automake to produce Makefile.in.
noinst_LTLIBRARIES = libmoleculizer_nauty.la

abs_top_builddir=@abs_top_builddir@
## SRCDIR = $(abs_top_builddir)/src
SRCDIR = ..

VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = moleculizer-$(VERSION)
libmoleculizer_nautydir=$(includedir)/$(MOLECULIZER_DIRNAME)/moleculizer/nauty

LIBFND = ../libmoleculizer/fnd/libmoleculizer_fnd.a
LIBUTL = ../libmoleculizer/utl/libmoleculizer_utl.a
LIBMZR = ../libmoleculizer/mzr/libmoleculizer_mzr.a
LIBSTOCH = ../libmoleculizer/stoch/libmoleculizer_stoch.a
LIBNMR = ../libmoleculizer/nmr/libmoleculizer_nmr.a
LIBCPX = ../libmoleculizer/cpx/libmoleculizer_cpx.a
LIBMOL = ../libmoleculizer/mol/libmoleculizer_mol.a
LIBPLEX = ../libmoleculizer/plex/libmoleculizer_plex.a
LIBDIMER = ../libmoleculizer/dimer/libmoleculizer_dimer.a
LIBFTR = ../libmoleculizer/ftr/libmoleculizer_ftr.a
LIBNAUTY = ../libmoleculizer/nauty/libmoleculizer_nauty.a

ALL_MZR_LIBS = $(LIBFND) $(LIBUTL) $(LIBMZR) $(LIBSTOCH) $(LIBNMR) $(LIBNMR) $(LIBCPX) $(LIBMOL) $(LIBPLEX) $(LIBDIMER) $(LIBFTR)

libmoleculizer_nauty_la_SOURCES =\
nauty.c \
nautil.c \
naugraph.c

EXTRA_DIST = gtools.h \
gutils.h \
naugroup.h \
nausparse.h \
nautaux.h \
nautinv.h \
naututil.h \
nauty.h \
planarity.h \
rng.h