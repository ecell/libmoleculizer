## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_utl.la

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

LIBXMLPP_CFLAGS = @LIBXMLPP_CFLAGS@
LIBXMLPP_LIBS = @LIBXMLPP_LIBS@

AM_CXXFLAGS = $(LIBXMLPP_CFLAGS) -I.. -I../../src @PYTHON_CPPFLAGS@

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
