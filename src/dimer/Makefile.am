## Process this file with automake to produce Makefile.in.

lib_LTLIBRARIES = libmoleculizer_dimer.la

abs_top_builddir=@abs_top_builddir@

SRCDIR = ..

VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = moleculizer-$(VERSION)
libmoleculizer_dimerdir=$(includedir)/$(MOLECULIZER_DIRNAME)/moleculizer/dimer

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

libmoleculizer_dimer_la_LIBADD = -lstdc++ $(LIBUTL) $(LIBFND) $(LIBCPX) $(LIBMZR) $(LIBMOL) $(LIBPLEX)
AM_CXXFLAGS = -I/usr/include/libxml++-1.0

libmoleculizer_dimer_la_SOURCES =\
decompRxnGen.cc \
decomposeExtrap.cc \
dimerDomParse.cc \
dimerEltName.cc \
dimerUnit.cc \
dimerizeExtrap.cc \
dimerizeRxnGen.cc \
missingDecomposeRateXcpt.cc \
missingDimerizeInvariantXcpt.cc \
missingDimerizeRateXcpt.cc

libmoleculizer_dimer_HEADERS=\
decompFam.hh \
decompRxnGen.hh \
decomposeExtrap.hh \
dimerEltName.hh \
dimerUnit.hh \
dimerizeExtrap.hh \
dimerizeFam.hh \
dimerizeRxnGen.hh \
missingDecomposeRateXcpt.hh \
missingDimerizeInvariantXcpt.hh \
missingDimerizeRateXcpt.hh 



