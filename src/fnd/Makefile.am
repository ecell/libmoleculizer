## Process this file with automake to produce Makefile.in.

lib_LTLIBRARIES = libmoleculizer_fnd.la

abs_top_builddir=@abs_top_builddir@

# SRCDIR = $(abs_top_builddir)/src
SRCDIR = ..

VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = moleculizer-$(VERSION)
libmoleculizer_fnddir=$(includedir)/$(MOLECULIZER_DIRNAME)/moleculizer/fnd

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


libmoleculizer_fnd_la_LIBADD = -lstdc++ 

AM_CXXFLAGS = -I/usr/include/libxml++-1.0

libmoleculizer_fnd_la_SOURCES =\
badDumpFileXcpt.cc \
dumpStream.cc \
pchem.cc \
physConst.cc \
speciesNotMassiveXcpt.cc

libmoleculizer_fnd_HEADERS=\
badDumpFileXcpt.hh \
basicDmpColumn.hh \
basicDumpable.hh \
basicReaction.hh \
basicSpecies.hh \
binaryRxnGen.hh \
dmpColumn.hh \
dumpStream.hh \
dumpable.hh \
event.hh \
feature.hh \
featureContext.hh \
featureMap.hh \
featureStimulus.hh \
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
reactionNetworkCatalog.hh \
reactionNetworkComponent.hh \
rxnGen.hh \
sensitive.hh \
sensitiveOnce.hh \
sensitivityList.hh \
speciesNotMassiveXcpt.hh \
stateVar.hh \
varDumpable.hh

