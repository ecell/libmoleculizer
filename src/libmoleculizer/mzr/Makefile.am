## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_mzr.la

VERSION = @INSTALLATION_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_mzrdir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/mzr

AM_CXXFLAGS = -Wall

libmoleculizer_mzr_la_CXXFLAGS = $(AM_CXXFLAGS) @LIBXMLPP_CFLAGS@ @PYTHON_CSPEC@
libmoleculizer_mzr_la_LDFLAGS = @PYTHON_LSPEC@

libmoleculizer_mzr_la_SOURCES =\
dumpUtils.cc \
libmzr_c_interface.cc \
moleculizer.cc \
mzrEltName.cc \
mzrException.cc \
mzrReaction.cc \
mzrSpecies.cc \
mzrUnit.cc \
mzrUnitInsert.cc \
mzrUnitParse.cc \
pythonRulesManager.cc \
spatialExtrapolationFunctions.cc \
unit.cc \
unitsMgr.cc

libmoleculizer_mzr_HEADERS=\
createEvent.hh \
dumpUtils.hh \
inputCapTest.hh \
libmzr_c_interface.h \
molarFactor.hh \
moleculizer.hh \
mzrEltName.hh \
mzrEvent.hh \
mzrException.hh \
mzrHelperFunctions.hh \
mzrReaction.hh \
mzrSpecies.hh \
mzrSpeciesDumpable.hh \
mzrSpeciesDumpableImpl.hh \
mzrStream.hh \
mzrUnit.hh \
pythonRulesManager.hh \
respondReaction.hh \
rxnDescriptionInterface.hh \
spatialExtrapolationFunctions.hh \
unit.hh \
unitsMgr.hh





