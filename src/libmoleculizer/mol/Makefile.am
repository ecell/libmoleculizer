## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_mol.la

abs_top_builddir=@abs_top_builddir@
abs_top_srcdir=@abs_top_srcdir@

## SRCDIR = $(abs_top_builddir)/src
SRCDIR = ..

VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_moldir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/mol

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

# libmoleculizer_mol_la_LIBADD = -lstdc++ $(LIBMZR) $(LIBCPX)

LIBXMLPP_CFLAGS = @LIBXMLPP_CFLAGS@
LIBXMLPP_LIBS = @LIBXMLPP_LIBS@

GSL_CFLAGS = @GSL_CFLAGS@
GSL_LIBS = @GSL_LIBS@

AM_CXXFLAGS = @CXXFLAGS@ -Wall $(LIBXMLPP_CFLAGS) @PYTHON_CPPFLAGS@ -I$(abs_top_srcdir)/src

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

