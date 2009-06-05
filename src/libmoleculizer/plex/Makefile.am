## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_plex.la

abs_top_builddir=@abs_top_builddir@
abs_top_srcdir=@abs_top_srcdir@

## SRCDIR = $(abs_top_builddir)/src
SRCDIR = ..

VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_plexdir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/plex

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

# libmoleculizer_plex_la_LIBADD = -lstdc++ $(LIBUTL) $(LIBNMR) $(LIBFND) $(LIBCPX) $(LIBMZR) $(LIBMOL) 

LIBXMLPP_CFLAGS = @LIBXMLPP_CFLAGS@
LIBXMLPP_LIBS = @LIBXMLPP_LIBS@

GSL_CFLAGS = @GSL_CFLAGS@
GSL_LIBS = @GSL_LIBS@

AM_CXXFLAGS = @CXXFLAGS@ -Wall $(LIBXMLPP_CFLAGS) @PYTHON_CPPFLAGS@ -I$(abs_top_srcdir)/src

if BUILD_LIBXMLPP
AM_CXXFLAGS += -I../../../libxml++-1.0.5
endif

libmoleculizer_plex_la_SOURCES =\
dupNodeOmniXcpt.cc \
multBoundSiteXcpt.cc \
mzrPlex.cc \
mzrPlexFamily.cc \
mzrPlexSpecies.cc \
mzrRecognizer.cc \
noKineticConstsXcpt.cc \
noOmniForNodeXcpt.cc \
parseOmniPlex.cc \
parsePlex.cc \
parsedPlexNotConnXcpt.cc \
parserPlex.cc \
plexEltName.cc \
plexUnit.cc \
plexUnitParse.cc \
unkMolInstXcpt.cc

libmoleculizer_plex_HEADERS=\
dupNodeOmniXcpt.hh \
multBoundSiteXcpt.hh \
mzrOmniPlex.hh \
mzrOmniStructureQueries.hh \
mzrPlex.hh \
mzrPlexFamily.hh \
mzrPlexQueries.hh \
mzrPlexSpecies.hh \
mzrRecognizer.hh \
noKineticConstsXcpt.hh \
noOmniForNodeXcpt.hh \
parseOmniPlex.hh \
parsePlex.hh \
parsedPlexNotConnXcpt.hh \
parserPlex.hh \
plexEltName.hh \
plexExceptions.hh \
plexUnit.hh \
unkMolInstXcpt.hh

