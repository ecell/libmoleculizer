## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_cpx.la

abs_top_builddir=@abs_top_builddir@

SRCDIR = ..

VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = moleculizer-$(VERSION)
libmoleculizer_cpxdir=$(includedir)/$(MOLECULIZER_DIRNAME)/moleculizer/cpx

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

LIBXML_CFLAGS = @LIBXML_CFLAGS@
LIBXML_LIBS = @LIBXML_LIBS@

GSL_CFLAGS = @GSL_CFLAGS@
GSL_LIBS = @GSL_LIBS@

AM_CXXFLAGS = @CXXFLAGS@ -Wall $(LIBXML_CFLAGS)

if STATIC_BUILD
AM_LDFLAGS= -all-static
endif

# libmoleculizer_cpx_la_LIBADD = $(LIBFND)

libmoleculizer_cpx_la_SOURCES =\
basicBndSite.cc \
cpxXcpt.cc \
modMolMixin.cc \
modStateMixin.cc \
plexIso.cc \
plexMap.cc \
siteToShapeMap.cc

libmoleculizer_cpx_HEADERS=\
alloMol.hh \
alloMolImpl.hh \
basicBndSite.hh \
basicMol.hh \
basicMolImpl.hh \
basicPlex.hh \
basicPlexImpl.hh \
binding.hh \
cpxXcpt.hh \
cxBinding.hh \
cxBindingImpl.hh \
cxMol.hh \
cxMolImpl.hh \
cxOmni.hh \
cxOmniImpl.hh \
cxSite.hh \
cxSiteImpl.hh \
ftrSpec.hh \
hashMolRec.hh \
hashMolRecImpl.hh \
isoSearch.hh \
isoSearchImpl.hh \
knownBindings.hh \
modMixinQuery.hh \
modMol.hh \
modMolImpl.hh \
modMolMixin.hh \
modMolState.hh \
modMolStateQuery.hh \
modStateMixin.hh \
modification.hh \
molState.hh \
molStateQuery.hh \
omniPlex.hh \
omniPlexFeature.hh \
omniStructureQuery.hh \
omniStructureQueryImpl.hh \
plexFamily.hh \
plexFamilyImpl.hh \
plexIso.hh \
plexIsoImpl.hh \
plexMap.hh \
plexMapImpl.hh \
plexQuery.hh \
plexQueryImpl.hh \
plexSpcsMixin.hh \
plexSpcsMixinCanonicalNamingImpl.hh \
plexSpcsMixinImpl.hh \
prm.hh \
queryAlloList.hh \
queryAlloListImpl.hh \
recognizer.hh \
recognizerImpl.hh \
reportIsoSearch.hh \
siteShape.hh \
siteToShapeMap.hh \
siteToShapeMapImpl.hh \
smallMol.hh \
stateMol.hh \
structuralBinding.hh \
subPlexSpec.hh

