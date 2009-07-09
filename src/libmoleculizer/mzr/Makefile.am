## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_mzr.la

abs_top_builddir=@abs_top_builddir@
abs_top_srcdir=@abs_top_srcdir@

SRCDIR = ..

VERSION = @PACKAGE_VERSION@

MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_mzrdir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/mzr

LIBXMLPP_CFLAGS = @LIBXMLPP_CFLAGS@
LIBXMLPP_LIBS = @LIBXMLPP_LIBS@

AM_CXXFLAGS = @CXXFLAGS@ -Wall $(LIBXMLPP_CFLAGS) -I.. @PYTHON_CSPEC@ -I$(abs_top_srcdir)/src

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





