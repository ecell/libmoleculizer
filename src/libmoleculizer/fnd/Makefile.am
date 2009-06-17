## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_fnd.la

abs_top_builddir=@abs_top_builddir@
abs_top_srcdir=@abs_top_srcdir@

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

VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_fnddir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/fnd

AM_CXXFLAGS = 

# libmoleculizer_fnd_la_LIBADD = -lstdc++ 

libmoleculizer_fnd_la_SOURCES =\
dumpStream.cc \
fndXcpt.cc \
pchem.cc \
physConst.cc

libmoleculizer_fnd_HEADERS=\
basicDmpColumn.hh \
basicDumpable.hh \
basicReaction.hh \
basicSpecies.hh \
binaryRxnGen.hh \
coreRxnGen.hh \
dmpColumn.hh \
dumpStream.hh \
dumpable.hh \
event.hh \
feature.hh \
featureContext.hh \
featureMap.hh \
featureStimulus.hh \
fndXcpt.hh \
gillspReaction.hh \
massive.hh \
multiSpeciesDumpable.hh \
newContextStimulus.hh \
newSpeciesStimulus.hh \
notifier.hh \
pchem.hh \
physConst.hh \
query.hh \
queryImpl.hh \
querySpeciesDumpable.hh \
reactionNetworkDescription.hh \
reactionNetworkDescriptionImpl.hh \
reactionNetworkComponent.hh \
rxnGen.hh \
sensitive.hh \
sensitiveOnce.hh \
sensitivityList.hh \
stateVar.hh \
varDumpable.hh

