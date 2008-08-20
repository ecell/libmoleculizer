## Process this file with automake to produce Makefile.in.
noinst_LIBRARIES = libmoleculizer_nauty.a

abs_top_builddir=@abs_top_builddir@
## SRCDIR = $(abs_top_builddir)/src
SRCDIR = ..

VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = moleculizer-$(VERSION)
libmoleculizer_nautydir=$(includedir)/$(MOLECULIZER_DIRNAME)/moleculizer/nauty

LIBFND = ../fnd/libmoleculizer_fnd.la
LIBUTL = ../utl/libmoleculizer_utl.la
LIBMZR = ../mzr/libmoleculizer_mzr.la
LIBSTOCH = ../stoch/libmoleculizer_stoch.la
LIBNMR = ../nmr/libmoleculizer_nmr.la
LIBCPX = ../cpx/libmoleculizer_cpx.la
LIBMOL = ../mol/libmoleculizer_mol.la
LIBPLEX = ../plex/libmoleculizer_plex.la
LIBDIMER = ../dimer/libmoleculizer_dimer.la
LIBFTR = ../ftr/libmoleculizer_ftr.la
LIBNAUTY = ../nauty/libmoleculizer_nauty.a
ALL_MZR_LIBS = $(LIBFND) $(LIBUTL) $(LIBMZR) $(LIBSTOCH) $(LIBNMR) $(LIBNMR) $(LIBCPX) $(LIBMOL) $(LIBPLEX) $(LIBDIMER) $(LIBFTR)

if STATIC_BUILD
AM_LDFLAGS= -all-static
endif

libmoleculizer_nauty_a_SOURCES =\
nauty.c \
nautil.c \
naugraph.c

libmoleculizer_nauty_HEADERS= \
gtools.h \
gutils.h \
naugroup.h \
nausparse.h \
nautaux.h \
nautinv.h \
naututil.h \
nauty.h	\
planarity.h \
rng.h