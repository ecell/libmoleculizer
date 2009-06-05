## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_ftr.la

abs_top_builddir=@abs_top_builddir@
abs_top_srcdir=@abs_top_srcdir@

SRCDIR = ..

VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_ftrdir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/ftr

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

MOLECULIZER_INCLUDE_DIR=${includedir}/moleculizer

LIBXMLPP_CFLAGS = @LIBXMLPP_CFLAGS@
LIBXMLPP_LIBS = @LIBXMLPP_LIBS@

AM_CXXFLAGS = @CXXFLAGS@ -Wall $(LIBXMLPP_CFLAGS) @PYTHON_CPPFLAGS@ -I$(abs_top_srcdir)/src

if BUILD_LIBXMLPP
AM_CXXFLAGS += -I../../../libxml++-1.0.5
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

