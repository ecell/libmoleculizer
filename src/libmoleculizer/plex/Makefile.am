## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_plex.la

VERSION = @INSTALLATION_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_plexdir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/plex

AM_CXXFLAGS = -Wall

libmoleculizer_plex_la_CXXFLAGS = $(AM_CXXFLAGS) @LIBXMLPP_CFLAGS@ @PYTHON_CSPEC@ 
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

