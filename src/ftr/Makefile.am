## Process this file with automake to produce Makefile.in.

lib_LTLIBRARIES = libmoleculizer_ftr.la

abs_top_builddir=@abs_top_builddir@

SRCDIR = ..

VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = moleculizer-$(VERSION)
libmoleculizer_ftrdir=$(includedir)/$(MOLECULIZER_DIRNAME)/moleculizer/ftr

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

MOLECULIZER_INCLUDE_DIR=${includedir}/moleculizer

libmoleculizer_ftr_la_LIBADD = -lstdc++ $(LIBUTL) $(LIBFND) $(LIBMZR) $(LIBMOL) $(LIBPLEX)

AM_CXXFLAGS = -I/usr/include/libxml++-1.0

if STATIC_BUILD
AM_LDFLAGS= -all-static
endif

libmoleculizer_ftr_la_SOURCES =\
badModMolInstanceXcpt.cc \
badSmallMolInstanceXcpt.cc \
badSmallMolXcpt.cc \
exchangedModXcpt.cc \
ftrDomParse.cc \
ftrEltName.cc \
omniExtrap.cc \
omniGen.cc \
parseOmniGen.cc \
parseUniMolGen.cc \
uniMolExtrap.cc \
uniMolGen.cc

libmoleculizer_ftr_HEADERS=\
badModMolInstanceXcpt.hh \
badSmallMolInstanceXcpt.hh \
badSmallMolXcpt.hh \
exchangedModXcpt.hh \
ftrEltName.hh \
ftrUnit.hh \
omniExtrap.hh \
omniFam.hh \
omniGen.hh \
parseOmniGen.hh \
parseUniMolGen.hh \
uniMolExtrap.hh \
uniMolFam.hh \
uniMolGen.hh

