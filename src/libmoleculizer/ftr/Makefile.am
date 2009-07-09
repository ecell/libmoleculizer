## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_ftr.la

VERSION = @INSTALLATION_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_ftrdir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/ftr

AM_CXXFLAGS = @LIBXMLPP_CFLAGS@ @PYTHON_CSPEC@ 

libmoleculizer_ftr_la_SOURCES =\
badModMolInstanceXcpt.cc \
badSmallMolInstanceXcpt.cc \
badSmallMolXcpt.cc \
exchangedModXcpt.cc \
ftrDomParse.cc \
ftrEltName.cc \
omniExtrap.cc \
omniGen.cc \
parseOmniGen.cc \
parseUniMolGen.cc \
uniMolExtrap.cc \
uniMolGen.cc

libmoleculizer_ftr_HEADERS=\
badModMolInstanceXcpt.hh \
badSmallMolInstanceXcpt.hh \
badSmallMolXcpt.hh \
exchangedModXcpt.hh \
ftrEltName.hh \
ftrUnit.hh \
omniExtrap.hh \
omniFam.hh \
omniGen.hh \
parseOmniGen.hh \
parseUniMolGen.hh \
uniMolExtrap.hh \
uniMolFam.hh \
uniMolGen.hh

