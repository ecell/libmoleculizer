## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_stoch.la

VERSION = @INSTALLATION_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_stochdir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/stoch

AM_CXXFLAGS = @LIBXMLPP_CFLAGS@ @PYTHON_CSPEC@ 

libmoleculizer_stoch_la_SOURCES =\
badStochSpeciesTagXcpt.cc \
stochDomParse.cc \
stochEltName.cc \
stochSpecies.cc \
stochUnit.cc \
unkStochSpeciesXcpt.cc

libmoleculizer_stoch_HEADERS=\
badStochSpeciesTagXcpt.hh \
stochEltName.hh \
stochSpecies.hh \
stochUnit.hh \
unkStochSpeciesXcpt.hh

