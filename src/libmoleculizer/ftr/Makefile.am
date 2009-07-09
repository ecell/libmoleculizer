## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_ftr.la

abs_top_builddir=@abs_top_builddir@
abs_top_srcdir=@abs_top_srcdir@

SRCDIR = ..

VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_ftrdir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/ftr

MOLECULIZER_INCLUDE_DIR=${includedir}/moleculizer

LIBXMLPP_CFLAGS = @LIBXMLPP_CFLAGS@
LIBXMLPP_LIBS = @LIBXMLPP_LIBS@

AM_CXXFLAGS = @CXXFLAGS@ -Wall $(LIBXMLPP_CFLAGS) @PYTHON_CSPEC@ -I$(abs_top_srcdir)/src

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

