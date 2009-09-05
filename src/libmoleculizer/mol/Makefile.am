## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_mol.la

VERSION = @INSTALLATION_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_moldir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/mol

AM_CXXFLAGS = -Wall

libmoleculizer_mol_la_CXXFLAGS = $(AM_CXXFLAGS) @LIBXMLPP_CFLAGS@ @PYTHON_CSPEC@ 
libmoleculizer_mol_la_SOURCES =\
badModMolXcpt.cpp \
badMolParamXcpt.cpp \
badSmallMolXcpt.cpp \
dupModNameXcpt.cpp \
dupModSiteNameXcpt.cpp \
dupMolNameXcpt.cpp \
dupSiteNameXcpt.cpp \
molDomInsert.cpp \
molDomParse.cpp \
molEltName.cpp \
molUnit.cpp \
mzrBndSite.cpp \
mzrModMol.cpp \
mzrMol.cpp \
smallMol.cpp \
unkModSiteXcpt.cpp \
unkModXcpt.cpp \
unkMolXcpt.cpp \
unkSiteShapeXcpt.cpp \
unkSiteXcpt.cpp

libmoleculizer_mol_HEADERS=\
badModMolXcpt.hpp \
badMolParamXcpt.hpp \
badSmallMolXcpt.hpp \
dupModNameXcpt.hpp \
dupModSiteNameXcpt.hpp \
dupMolNameXcpt.hpp \
dupSiteNameXcpt.hpp \
molEltName.hpp \
molUnit.hpp \
molXcpt.hpp \
mzrBndSite.hpp \
mzrModMol.hpp \
mzrMol.hpp \
parseBndSite.hpp \
parseMod.hpp \
parseModMap.hpp \
parseModSite.hpp \
parseSiteShapeName.hpp \
siteFeature.hpp \
smallMol.hpp \
unkModSiteXcpt.hpp \
unkModXcpt.hpp \
unkMolXcpt.hpp \
unkSiteShapeXcpt.hpp \
unkSiteXcpt.hpp

