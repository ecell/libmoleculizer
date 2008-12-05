## Process this file with automake to produce Makefile.in.

lib_LTLIBRARIES = libmoleculizer_fnd.la

abs_top_builddir=@abs_top_builddir@

## SRCDIR = $(abs_top_builddir)/src
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


LIBXML_CFLAGS = @LIBXML_CFLAGS@
LIBXML_LIBS = @LIBXML_LIBS@

GSL_CFLAGS = @GSL_CFLAGS@
GSL_LIBS = @GSL_LIBS@

AM_CXXFLAGS = @CXXFLAGS@ -Wall $(LIBXML_CFLAGS)

if STATIC_BUILD
AM_LDFLAGS= -all-static
endif


libmoleculizer_fnd_la_LIBADD = -lstdc++ 

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
reactionNetworkComponent.hh \
rxnGen.hh \
sensitive.hh \
sensitiveOnce.hh \
sensitivityList.hh \
stateVar.hh \
utility.hh \
varDumpable.hh

