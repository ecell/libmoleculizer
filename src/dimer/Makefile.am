## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_dimer.la

abs_top_builddir=@abs_top_builddir@
abs_top_srcdir=@abs_top_srcdir@

SRCDIR = ..

VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_dimerdir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/dimer

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

# libmoleculizer_dimer_la_LIBADD = -lstdc++ $(LIBUTL) $(LIBFND) $(LIBCPX) $(LIBMZR) $(LIBMOL) $(LIBPLEX)

LIBXMLPP_CFLAGS = @LIBXMLPP_CFLAGS@
LIBXMLPP_LIBS = @LIBXMLPP_LIBS@

GSL_CFLAGS = @GSL_CFLAGS@
GSL_LIBS = @GSL_LIBS@

AM_CXXFLAGS = @CXXFLAGS@ -Wall $(LIBXMLPP_CFLAGS) @PYTHON_CPPFLAGS@ -I$(abs_top_srcdir)/src

libmoleculizer_dimer_la_SOURCES =\
decompRxnGen.cc \
decomposeExtrap.cc \
dimerDomParse.cc \
dimerEltName.cc \
dimerUnit.cc \
dimerXcpt.cc \
dimerizeExtrap.cc \
dimerizeRxnGen.cc


libmoleculizer_dimer_HEADERS=\
decompFam.hh \
decompRxnGen.hh \
decomposeExtrap.hh \
dimerEltName.hh \
dimerUnit.hh \
dimerXcpt.hh \
dimerizeExtrap.hh \
dimerizeFam.hh \
dimerizeRxnGen.hh



