## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_utl.la

VERSION = @INSTALLATION_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_utldir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/utl

AM_CXXFLAGS = -Wall

libmoleculizer_utl_la_CXXFLAGS = $(AM_CXXFLAGS) @LIBXMLPP_CFLAGS@
libmoleculizer_utl_la_LDFLAGS = @LIBXMLPP_LIBS@
libmoleculizer_utl_la_SOURCES =\
arg.cpp \
dom.cpp \
domXcpt.cpp \
frexp10.cpp \
linearHash.cpp \
utlXcpt.cpp \
utility.cpp \
utlEltName.cpp

libmoleculizer_utl_HEADERS=\
arg.hpp \
AssocVector.h \
autoCache.hpp \
autoCatalog.hpp \
autoVector.hpp \
defs.hpp \
dom.hpp \
domJob.hpp \
domJobImpl.hpp \
domXcpt.hpp \
forBoth.hpp \
forceInsert.hpp \
frexp10.hpp \
funcInsert.hpp \
linearHash.hpp \
message.hpp \
utility.hpp \
utlEltName.hpp \
utlHelper.hpp \
utlXcpt.hpp \
warn.hpp \
writeOutputGraph.hpp \
xcpt.hpp
