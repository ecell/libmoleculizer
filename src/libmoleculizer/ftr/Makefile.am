## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_ftr.la

VERSION = @INSTALLATION_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_ftrdir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/ftr

AM_CXXFLAGS = -Wall

libmoleculizer_ftr_la_CXXFLAGS = $(AM_CXXFLAGS) @LIBXMLPP_CFLAGS@ @PYTHON_CSPEC@ 
libmoleculizer_ftr_la_SOURCES =\
badModMolInstanceXcpt.cpp \
badSmallMolInstanceXcpt.cpp \
badSmallMolXcpt.cpp \
exchangedModXcpt.cpp \
ftrDomParse.cpp \
ftrEltName.cpp \
omniExtrap.cpp \
omniGen.cpp \
parseOmniGen.cpp \
parseUniMolGen.cpp \
uniMolExtrap.cpp \
uniMolGen.cpp

libmoleculizer_ftr_HEADERS=\
badModMolInstanceXcpt.hpp \
badSmallMolInstanceXcpt.hpp \
badSmallMolXcpt.hpp \
exchangedModXcpt.hpp \
ftrEltName.hpp \
ftrUnit.hpp \
omniExtrap.hpp \
omniFam.hpp \
omniGen.hpp \
parseOmniGen.hpp \
parseUniMolGen.hpp \
uniMolExtrap.hpp \
uniMolFam.hpp \
uniMolGen.hpp

