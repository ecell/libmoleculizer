## Process this file with automake to produce Makefile.in.

lib_LTLIBRARIES = libmoleculizer_mzr.la

abs_top_builddir=@abs_top_builddir@

## SRCDIR = $(abs_top_builddir)/src
SRCDIR = ..


VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = moleculizer-$(VERSION)
libmoleculizer_mzrdir=$(includedir)/$(MOLECULIZER_DIRNAME)/moleculizer/mzr

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

AM_CXXFLAGS = -I/usr/include/libxml++-1.0

libmoleculizer_mzr_la_LIBADD = -lstdc++ $(LIBUTL)

libmoleculizer_mzr_la_SOURCES =\
continuator.cc \
createEvent.cc \
dumpStateEvent.cc \
dumpUtils.cc \
dumpableNotSpeciesStreamXcpt.cc \
dupDumpableNameXcpt.cc \
dupSpeciesNameXcpt.cc \
growEvent.cc \
moleculizer.cc \
moleculizerDebugFunctions.cc \
mzrEltName.cc \
mzrReaction.cc \
mzrUnit.cc \
mzrUnitInsert.cc \
mzrUnitParse.cc \
noReactEvent.cc \
noStopEventWarning.cc \
species.cc \
stopEvent.cc \
stopEventInPastXcpt.cc \
tabDumpEvent.cc \
unitsMgr.cc \
unkDumpableXcpt.cc \
unkSpeciesXcpt.cc \
unkStatStreamXcpt.cc \
volumeEvent.cc

libmoleculizer_mzr_HEADERS=\
debug.hh \
continuator.hh \
createEvent.hh \
dumpStateEvent.hh \
dumpUtils.hh \
dumpableNotSpeciesStreamXcpt.hh \
dupDumpableNameXcpt.hh \
dupSpeciesNameXcpt.hh \
eventQueue.hh \
growEvent.hh \
inputCapTest.hh \
inputCapXcpt.hh \
molarFactor.hh \
moleculizer.hh \
mzrEltName.hh \
mzrEvent.hh \
mzrReaction.hh \
mzrSpecies.hh \
mzrSpeciesDumpable.hh \
mzrSpeciesDumpableImpl.hh \
mzrStream.hh \
mzrUnit.hh \
noReactEvent.hh \
noStopEventWarning.hh \
reactCountDumpable.hh \
respondReaction.hh \
rxnDescriptionInterface.hh \
simTimeDumpable.hh \
stopEvent.hh \
stopEventInPastXcpt.hh \
tabDumpEvent.hh \
unit.hh \
unitsMgr.hh \
unkDumpableXcpt.hh \
unkSpeciesXcpt.hh \
unkStatStreamXcpt.hh \
volumeEvent.hh

