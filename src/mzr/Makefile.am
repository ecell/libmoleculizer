## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_mzr.la

abs_top_builddir=@abs_top_builddir@
abs_top_srcdir=@abs_top_srcdir@

SRCDIR = ..
SUBDIRS = tests


VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_mzrdir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/mzr

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

AM_CXXFLAGS = @CXXFLAGS@ -Wall $(LIBXMLPP_CFLAGS) -I.. @PYTHON_CPPFLAGS@ -I$(abs_top_srcdir)/src

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





