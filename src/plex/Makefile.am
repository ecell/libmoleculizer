lib_LTLIBRARIES = libmoleculizer_plex.la

abs_top_builddir=@abs_top_builddir@

# SRCDIR = $(abs_top_builddir)/src
SRCDIR = ..

VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = moleculizer-$(VERSION)
libmoleculizer_plexdir=$(includedir)/$(MOLECULIZER_DIRNAME)/moleculizer/plex

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

libmoleculizer_plex_la_LIBADD = -lstdc++ $(LIBUTL) $(LIBNMR) $(LIBFND) $(LIBCPX) $(LIBMZR) $(LIBMOL)
AM_CXXFLAGS = -I/usr/include/libxml++-1.0

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
plexUnit.hh \
unkMolInstXcpt.hh

