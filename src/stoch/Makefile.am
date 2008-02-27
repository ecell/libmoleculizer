## Process this file with automake to produce Makefile.in.

lib_LTLIBRARIES = libmoleculizer_stoch.la

abs_top_builddir=@abs_top_builddir@

## SRCDIR = $(abs_top_builddir)/src
SRCDIR = ..

VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = moleculizer-$(VERSION)
libmoleculizer_stochdir=$(includedir)/$(MOLECULIZER_DIRNAME)/moleculizer/stoch

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
ALL_MZR_LIBS = $(LIBFND) $(LIBUTL) $(LIBMZR) $(LIBSTOCH) $(LIBNMR) $(LIBNMR) $(LIBCPX) $(LIBMOL) $(LIBPLEX) $(LIBDIMER) $(LIBFTR)


libmoleculizer_stoch_la_LIBADD = -lstdc++ $(LIBUTL) $(LIBFND) $(LIBMZR)
AM_CXXFLAGS = -I/usr/include/libxml++-1.0

if STATIC_BUILD
AM_CXXFLAGS += -pg
AM_LDFLAGS= -all-static
endif

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

