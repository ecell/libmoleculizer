## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_mzr.la

VERSION = @INSTALLATION_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_mzrdir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/mzr

AM_CXXFLAGS = -Wall

libmoleculizer_mzr_la_CXXFLAGS = $(AM_CXXFLAGS) @LIBXMLPP_CFLAGS@ @PYTHON_CSPEC@
libmoleculizer_mzr_la_LDFLAGS = @PYTHON_LSPEC@

libmoleculizer_mzr_la_SOURCES =\
dumpUtils.cpp \
moleculizer.cpp \
mzrEltName.cpp \
mzrException.cpp \
mzrReaction.cpp \
mzrSpecies.cpp \
mzrUnit.cpp \
mzrUnitInsert.cpp \
mzrUnitParse.cpp \
pythonRulesManager.cpp \
spatialExtrapolationFunctions.cpp \
unit.cpp \
unitsMgr.cpp

libmoleculizer_mzr_HEADERS=\
createEvent.hpp \
dumpUtils.hpp \
inputCapTest.hpp \
molarFactor.hpp \
moleculizer.hpp \
mzrEltName.hpp \
mzrEvent.hpp \
mzrException.hpp \
mzrHelperFunctions.hpp \
mzrReaction.hpp \
mzrSpecies.hpp \
mzrSpeciesDumpable.hpp \
mzrSpeciesDumpableImpl.hpp \
mzrStream.hpp \
mzrUnit.hpp \
pythonRulesManager.hpp \
respondReaction.hpp \
rxnDescriptionInterface.hpp \
spatialExtrapolationFunctions.hpp \
unit.hpp \
unitsMgr.hpp





