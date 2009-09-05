## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_stoch.la

VERSION = @INSTALLATION_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_stochdir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/stoch

AM_CXXFLAGS = -Wall

libmoleculizer_stoch_la_CXXFLAGS = $(AM_CXXFLAGS) @LIBXMLPP_CFLAGS@ @PYTHON_CSPEC@ 
libmoleculizer_stoch_la_SOURCES =\
badStochSpeciesTagXcpt.cpp \
stochDomParse.cpp \
stochEltName.cpp \
stochSpecies.cpp \
stochUnit.cpp \
unkStochSpeciesXcpt.cpp

libmoleculizer_stoch_HEADERS=\
badStochSpeciesTagXcpt.hpp \
stochEltName.hpp \
stochSpecies.hpp \
stochUnit.hpp \
unkStochSpeciesXcpt.hpp

