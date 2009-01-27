## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_utl.la

abs_top_builddir=@abs_top_builddir@

## SRCDIR = $(abs_top_builddir)/src
SRCDIR = ..

VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_utldir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/utl

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


if STATIC_BUILD
AM_LDFLAGS= -all-static
endif

LIBXML_CFLAGS = @LIBXML_CFLAGS@
LIBXML_LIBS = @LIBXML_LIBS@

AM_CXXFLAGS = $(LIBXML_CFLAGS) 

# libmoleculizer_utl_la_LIBADD = -lstdc++ $(LIBXML_LIBS) 

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
autoCache.hh \
autoCatalog.hh \
autoVector.hh \
defs.hh \
dom.hh \
domJob.hh \
domJobImpl.hh \
forBoth.hh \
forceInsert.hh \
frexp10.hh \
funcInsert.hh \
linearHash.hh \
message.hh \
utility.hh \
utlEltName.hh \
warn.hh \
xcpt.hh
