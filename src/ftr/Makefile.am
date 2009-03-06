## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_ftr.la

abs_top_builddir=@abs_top_builddir@

SRCDIR = ..

VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_ftrdir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/ftr

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

# libmoleculizer_ftr_la_LIBADD = -lstdc++ $(LIBUTL) $(LIBFND) $(LIBMZR) $(LIBMOL) $(LIBPLEX)

LIBXMLPP_CFLAGS = @LIBXMLPP_CFLAGS@
LIBXMLPP_LIBS = @LIBXMLPP_LIBS@

AM_CXXFLAGS = @CXXFLAGS@ -Wall $(LIBXMLPP_CFLAGS) @PYTHON_CPPFLAGS@

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

