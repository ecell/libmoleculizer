## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_fnd.la

abs_top_builddir=@abs_top_builddir@
abs_top_srcdir=@abs_top_srcdir@

## SRCDIR = $(abs_top_builddir)/src
SRCDIR = ..

VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_fnddir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/fnd

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


LIBXMLPP_CFLAGS = @LIBXMLPP_CFLAGS@
LIBXMLPP_LIBS = @LIBXMLPP_LIBS@

GSL_CFLAGS = @GSL_CFLAGS@
GSL_LIBS = @GSL_LIBS@

AM_CXXFLAGS = @CXXFLAGS@ -Wall $(LIBXMLPP_CFLAGS) -I$(abs_top_srcdir)/src

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

