## Process this file with automake to produce Makefile.in.

abs_top_srcdir=@abs_top_srcdir@
noinst_LTLIBRARIES = libmoleculizer_utl.la

LIBFND = $(abs_top_srcdir)/src/libmoleculizer/fnd/libmoleculizer_fnd.la
LIBUTL = $(abs_top_srcdir)/src/libmoleculizer/utl/libmoleculizer_utl.la
LIBMZR = $(abs_top_srcdir)/src/libmoleculizer/mzr/libmoleculizer_mzr.la
LIBSTOCH = $(abs_top_srcdir)/src/libmoleculizer/stoch/libmoleculizer_stoch.la
LIBNMR = $(abs_top_srcdir)/src/libmoleculizer/nmr/libmoleculizer_nmr.la
LIBCPX = $(abs_top_srcdir)/src/libmoleculizer/cpx/libmoleculizer_cpx.la
LIBMOL = $(abs_top_srcdir)/src/libmoleculizer/mol/libmoleculizer_mol.la
LIBPLEX = $(abs_top_srcdir)/src/libmoleculizer/plex/libmoleculizer_plex.la
LIBDIMER = $(abs_top_srcdir)/src/libmoleculizer/dimer/libmoleculizer_dimer.la
LIBFTR = $(abs_top_srcdir)/src/libmoleculizer/ftr/libmoleculizer_ftr.la
LIBNAUTY = $(abs_top_srcdir)/src/libmoleculizer/nauty/libmoleculizer_nauty.la


ALL_MZR_LIBS = $(LIBFND) $(LIBUTL) $(LIBMZR) $(LIBSTOCH) $(LIBNMR) $(LIBNMR) $(LIBCPX) $(LIBMOL) $(LIBPLEX) $(LIBDIMER) $(LIBFTR)

VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)

libmoleculizer_utldir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/utl

LIBXMLPP_CFLAGS = @LIBXMLPP_CFLAGS@
LIBXMLPP_LIBS = @LIBXMLPP_LIBS@

AM_CXXFLAGS = $(LIBXMLPP_CFLAGS)

if BUILD_LIBXMLPP
AM_CXXFLAGS += -I../../../libxml++-1.0.5
libmoleculizer_utl_la_LIBADD=../../../libxml++-1.0.5/libxml++/libxmlpp.la
endif

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
