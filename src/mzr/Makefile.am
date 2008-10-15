## Process this file with automake to produce Makefile.in.

lib_LTLIBRARIES = libmoleculizer_mzr.la

abs_top_builddir=@abs_top_builddir@

SRCDIR = ..
SUBDIRS = tests


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

LIBXML_CFLAGS = @LIBXML_CFLAGS@
LIBXML_LIBS = @LIBXML_LIBS@

AM_CXXFLAGS = @CXXFLAGS@ -Wall $(LIBXML_CFLAGS)

if STATIC_BUILD
AM_LDFLAGS= -all-static
endif

libmoleculizer_mzr_la_LIBADD = -lstdc++ $(LIBUTL)

libmoleculizer_mzr_la_SOURCES =\
dumpUtils.cc \
extrapolationFormula.cc \
moleculizer.cc \
mzrEltName.cc \
mzrException.cc \
mzrReaction.cc \
mzrSpecies.cc \
mzrUnit.cc \
mzrUnitInsert.cc \
mzrUnitParse.cc \
unitsMgr.cc

# growEvent.cc 

libmoleculizer_mzr_HEADERS=\
createEvent.hh \
dumpUtils.hh \
inputCapTest.hh \
molarFactor.hh \
moleculizer.hh \
mzrEltName.hh \
mzrEvent.hh \
mzrException.hh \
mzrHelperFunctions.hh \
mzrReaction.hh \
mzrSpecies.hh \
mzrStream.hh \
mzrUnit.hh \
respondReaction.hh \
rxnDescriptionInterface.hh \
unit.hh \
unitsMgr.hh





