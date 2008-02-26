## Process this file with automake to produce Makefile.in.

lib_LTLIBRARIES = libmoleculizer_utl.la

abs_top_builddir=@abs_top_builddir@

## SRCDIR = $(abs_top_builddir)/src
SRCDIR = ..

VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = moleculizer-$(VERSION)
libmoleculizer_utldir=$(includedir)/$(MOLECULIZER_DIRNAME)/moleculizer/utl

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
libmoleculizer_utl_la_LIBADD = -lstdc++ -lxml2 -lxml++-1.0 -lgsl -lgslcblas

libmoleculizer_utl_la_SOURCES =\
arg.cc \
badChildCountXcpt.cc \
badDoubleAttrXcpt.cc \
badElementCastXcpt.cc \
badIntAttrXcpt.cc \
badNNDoubleArgXcpt.cc \
badNNDoubleAttrXcpt.cc \
badNNIntArgXcpt.cc \
badNNIntAttrXcpt.cc \
badPosDoubleAttrXcpt.cc \
badPosIntArgXcpt.cc \
badPosIntAttrXcpt.cc \
dom.cc \
frexp10.cc \
gsl.cc \
insuffArgsXcpt.cc \
linearHash.cc \
missingAttrXcpt.cc \
string.cc \
unkArgXcpt.cc \
utility.cc \
utlEltName.cc

libmoleculizer_utl_HEADERS=\
arg.hh \
autoCache.hh \
autoCatalog.hh \
autoVector.hh \
badChildCountXcpt.hh \
badDoubleAttrXcpt.hh \
badElementCastXcpt.hh \
badFileNameXcpt.hh \
badIntAttrXcpt.hh \
badNNDoubleArgXcpt.hh \
badNNDoubleAttrXcpt.hh \
badNNIntArgXcpt.hh \
badNNIntAttrXcpt.hh \
badPosDoubleAttrXcpt.hh \
badPosIntArgXcpt.hh \
badPosIntAttrXcpt.hh \
debug.hh \
dom.hh \
domJob.hh \
domJobImpl.hh \
forBoth.hh \
forceInsert.hh \
frexp10.hh \
funcInsert.hh \
gsl.hh \
insuffArgsXcpt.hh \
linearHash.hh \
message.hh \
missingAttrXcpt.hh \
modelPreviouslyLoadedXcpt.hh \
noDocumentParsedXcpt.hh \
platform.hh \
string.hh \
unkArgXcpt.hh \
utility.hh \
utlEltName.hh \
warn.hh \
xcpt.hh
