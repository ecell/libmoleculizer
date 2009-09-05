## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_dimer.la

VERSION = @INSTALLATION_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_dimerdir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/dimer

AM_CXXFLAGS = -Wall 

libmoleculizer_dimer_la_CXXFLAGS = $(AM_CXXFLAGS) @LIBXMLPP_CFLAGS@ @PYTHON_CSPEC@ 
libmoleculizer_dimer_la_SOURCES =\
decompRxnGen.cpp \
decomposeExtrap.cpp \
dimerDomParse.cpp \
dimerEltName.cpp \
dimerUnit.cpp \
dimerXcpt.cpp \
dimerizeExtrap.cpp \
dimerizeRxnGen.cpp

libmoleculizer_dimer_HEADERS=\
decompFam.hpp \
decompRxnGen.hpp \
decomposeExtrap.hpp \
dimerEltName.hpp \
dimerUnit.hpp \
dimerXcpt.hpp \
dimerizeExtrap.hpp \
dimerizeFam.hpp \
dimerizeRxnGen.hpp



