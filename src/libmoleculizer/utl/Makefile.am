## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_utl.la

VERSION = @INSTALLATION_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_utldir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/utl

AM_CXXFLAGS = -Wall

libmoleculizer_utl_la_CXXFLAGS = $(AM_CXXFLAGS) @LIBXMLPP_CFLAGS@
libmoleculizer_utl_la_LDFLAGS = @LIBXMLPP_LIBS@
libmoleculizer_utl_la_SOURCES =\
arg.cc \
dom.cc \
domXcpt.cc \
frexp10.cc \
linearHash.cc \
utlXcpt.cc \
utility.cc \
utlEltName.cc

libmoleculizer_utl_HEADERS=\
arg.hh \
AssocVector.h \
autoCache.hh \
autoCatalog.hh \
autoVector.hh \
defs.hh \
dom.hh \
domJob.hh \
domJobImpl.hh \
domXcpt.hh \
forBoth.hh \
forceInsert.hh \
frexp10.hh \
funcInsert.hh \
linearHash.hh \
message.hh \
utility.hh \
utlEltName.hh \
utlHelper.hh \
utlXcpt.hh \
warn.hh \
writeOutputGraph.hh \
xcpt.hh
