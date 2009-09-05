## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_plex.la

VERSION = @INSTALLATION_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_plexdir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/plex

AM_CXXFLAGS = -Wall

libmoleculizer_plex_la_CXXFLAGS = $(AM_CXXFLAGS) @LIBXMLPP_CFLAGS@ @PYTHON_CSPEC@ 
libmoleculizer_plex_la_SOURCES =\
dupNodeOmniXcpt.cpp \
multBoundSiteXcpt.cpp \
mzrPlex.cpp \
mzrPlexFamily.cpp \
mzrPlexSpecies.cpp \
mzrRecognizer.cpp \
noKineticConstsXcpt.cpp \
noOmniForNodeXcpt.cpp \
parseOmniPlex.cpp \
parsePlex.cpp \
parsedPlexNotConnXcpt.cpp \
parserPlex.cpp \
plexEltName.cpp \
plexUnit.cpp \
plexUnitParse.cpp \
unkMolInstXcpt.cpp

libmoleculizer_plex_HEADERS=\
dupNodeOmniXcpt.hpp \
multBoundSiteXcpt.hpp \
mzrOmniPlex.hpp \
mzrOmniStructureQueries.hpp \
mzrPlex.hpp \
mzrPlexFamily.hpp \
mzrPlexQueries.hpp \
mzrPlexSpecies.hpp \
mzrRecognizer.hpp \
noKineticConstsXcpt.hpp \
noOmniForNodeXcpt.hpp \
parseOmniPlex.hpp \
parsePlex.hpp \
parsedPlexNotConnXcpt.hpp \
parserPlex.hpp \
plexEltName.hpp \
plexExceptions.hpp \
plexUnit.hpp \
unkMolInstXcpt.hpp

