## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_stoch.la

abs_top_builddir=@abs_top_builddir@
abs_top_srcdir=@abs_top_srcdir@

## SRCDIR = $(abs_top_builddir)/src
SRCDIR = ..

VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_stochdir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/stoch

LIBFND = $(abs_top_srcdir)/src/libmoleculizer/fnd/libmoleculizer_fnd.la
LIBUTL = $(abs_top_srcdir)/src/libmoleculizer/utl/libmoleculizer_utl.la
LIBMZR = $(abs_top_srcdir)/src/libmoleculizer/mzr/libmoleculizer_mzr.la
LIBSTOCH = $(abs_top_srcdir)/src/libmoleculizer/stoch/libmoleculizer_stoch.la
LIBNMR = $(abs_top_srcdir)/src/libmoleculizer/nmr/libmoleculizer_nmr.la
LIBCPX = $(abs_top_srcdir)/src/libmoleculizer/cpx/libmoleculizer_cpx.la
LIBMOL = $(abs_top_srcdir)/src/libmoleculizer/mol/libmoleculizer_mol.la
LIBPLEX = $(abs_top_srcdir)/src/libmoleculizer/plex/libmoleculizer_plex.la
LIBDIMER = $(abs_top_srcdir)/src/libmoleculizer/dimer/libmoleculizer_dimer.la
LIBFTR = $(abs_top_srcdir)/src/libmoleculizer/ftr/libmoleculizer_ftr.la
LIBNAUTY = $(abs_top_srcdir)/src/libmoleculizer/nauty/libmoleculizer_nauty.la


ALL_MZR_LIBS = $(LIBFND) $(LIBUTL) $(LIBMZR) $(LIBSTOCH) $(LIBNMR) $(LIBNMR) $(LIBCPX) $(LIBMOL) $(LIBPLEX) $(LIBDIMER) $(LIBFTR)


LIBXMLPP_CFLAGS = @LIBXMLPP_CFLAGS@
LIBXMLPP_LIBS = @LIBXMLPP_LIBS@

GSL_CFLAGS = @GSL_CFLAGS@
GSL_LIBS = @GSL_LIBS@

AM_CXXFLAGS = @CXXFLAGS@ -Wall $(LIBXMLPP_CFLAGS) @PYTHON_CSPEC@ -I$(abs_top_srcdir)/src

# libmoleculizer_stoch_la_LIBADD = -lstdc++ $(LIBUTL) $(LIBFND) $(LIBMZR)

libmoleculizer_stoch_la_SOURCES =\
badStochSpeciesTagXcpt.cc \
stochDomParse.cc \
stochEltName.cc \
stochSpecies.cc \
stochUnit.cc \
unkStochSpeciesXcpt.cc

libmoleculizer_stoch_HEADERS=\
badStochSpeciesTagXcpt.hh \
stochEltName.hh \
stochSpecies.hh \
stochUnit.hh \
unkStochSpeciesXcpt.hh

