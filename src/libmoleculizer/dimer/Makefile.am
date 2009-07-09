## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_dimer.la

VERSION = @INSTALLATION_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_dimerdir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/dimer

LIBXMLPP_CFLAGS = @LIBXMLPP_CFLAGS@
LIBXMLPP_LIBS = @LIBXMLPP_LIBS@

AM_CXXFLAGS = -Wall @LIBXMLPP_CFLAGS@ @PYTHON_CSPEC@ 

libmoleculizer_dimer_la_SOURCES =\
decompRxnGen.cc \
decomposeExtrap.cc \
dimerDomParse.cc \
dimerEltName.cc \
dimerUnit.cc \
dimerXcpt.cc \
dimerizeExtrap.cc \
dimerizeRxnGen.cc

libmoleculizer_dimer_HEADERS=\
decompFam.hh \
decompRxnGen.hh \
decomposeExtrap.hh \
dimerEltName.hh \
dimerUnit.hh \
dimerXcpt.hh \
dimerizeExtrap.hh \
dimerizeFam.hh \
dimerizeRxnGen.hh



