## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_mol.la

abs_top_builddir=@abs_top_builddir@
abs_top_srcdir=@abs_top_srcdir@

## SRCDIR = $(abs_top_builddir)/src
SRCDIR = ..

VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_moldir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/mol

LIBXMLPP_CFLAGS = @LIBXMLPP_CFLAGS@
LIBXMLPP_LIBS = @LIBXMLPP_LIBS@

GSL_CFLAGS = @GSL_CFLAGS@
GSL_LIBS = @GSL_LIBS@

AM_CXXFLAGS = @CXXFLAGS@ -Wall $(LIBXMLPP_CFLAGS) @PYTHON_CSPEC@ -I$(abs_top_srcdir)/src

libmoleculizer_mol_la_SOURCES =\
badModMolXcpt.cc \
badMolParamXcpt.cc \
badSmallMolXcpt.cc \
dupModNameXcpt.cc \
dupModSiteNameXcpt.cc \
dupMolNameXcpt.cc \
dupSiteNameXcpt.cc \
molDomInsert.cc \
molDomParse.cc \
molEltName.cc \
molUnit.cc \
mzrBndSite.cc \
mzrModMol.cc \
mzrMol.cc \
smallMol.cc \
unkModSiteXcpt.cc \
unkModXcpt.cc \
unkMolXcpt.cc \
unkSiteShapeXcpt.cc \
unkSiteXcpt.cc

libmoleculizer_mol_HEADERS=\
badModMolXcpt.hh \
badMolParamXcpt.hh \
badSmallMolXcpt.hh \
dupModNameXcpt.hh \
dupModSiteNameXcpt.hh \
dupMolNameXcpt.hh \
dupSiteNameXcpt.hh \
molEltName.hh \
molUnit.hh \
molXcpt.hh \
mzrBndSite.hh \
mzrModMol.hh \
mzrMol.hh \
parseBndSite.hh \
parseMod.hh \
parseModMap.hh \
parseModSite.hh \
parseSiteShapeName.hh \
siteFeature.hh \
smallMol.hh \
unkModSiteXcpt.hh \
unkModXcpt.hh \
unkMolXcpt.hh \
unkSiteShapeXcpt.hh \
unkSiteXcpt.hh

