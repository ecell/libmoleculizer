## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_cpx.la

abs_top_builddir=@abs_top_builddir@
abs_top_srcdir=@abs_top_srcdir@

SRCDIR = ..

VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_cpxdir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/cpx

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

LIBXMLPP_CFLAGS = @LIBXMLPP_CFLAGS@
LIBXMLPP_LIBS = @LIBXMLPP_LIBS@

GSL_CFLAGS = @GSL_CFLAGS@
GSL_LIBS = @GSL_LIBS@

AM_CXXFLAGS = @CXXFLAGS@ -Wall $(LIBXMLPP_CFLAGS) -I$(abs_top_srcdir)/src

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

