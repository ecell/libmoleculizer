## Process this file with automake to produce Makefile.in.

lib_LTLIBRARIES = libmoleculizer_nmr.la

abs_top_builddir=@abs_top_builddir@

SRCDIR = .. 

SUBDIRS = tests

VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = moleculizer-$(VERSION)
libmoleculizer_nmrdir=$(includedir)/$(MOLECULIZER_DIRNAME)/moleculizer/nmr

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

LIBXML_CFLAGS = @LIBXML_CFLAGS@
LIBXML_LIBS = @LIBXML_LIBS@

GSL_CFLAGS = @GSL_CFLAGS@
GSL_LIBS = @GSL_LIBS@

AM_CXXFLAGS = $(LIBXML_CFLAGS)

if STATIC_BUILD
AM_LDFLAGS = -all-static
endif

libmoleculizer_nmr_la_LIBADD = nauty/nautil.o nauty/nautinv.o nauty/nauty.o nauty/gtools.o nauty/naugraph.o

libmoleculizer_nmr_la_SOURCES =\
basicNameAssembler.cc \
complexOutputState.cc \
complexSpecies.cc \
complexSpeciesEncoderNames.cc \
complexSpeciesOutputMinimizer.cc \
mangledNameAssembler.cc \
namedMolecule.cc \
nmrEltName.cc \
nmrUnit.cc \
partialTokenList.cc \
permutation.cc \
permutationName.cc \
readableNameAssembler.cc 

libmoleculizer_nmr_HEADERS=\
abstractMol.hh \
basicNameAssembler.hh \
complexOutputState.hh \
complexSpeciesEncoderNames.hh \
complexSpecies.hh \
complexSpeciesOutputMinimizer.hh \
mangledNameAssembler.hh \
nameAssembler.hh \
nameAssemblers.hh \
namedMolecule.hh \
nameEncoderFactory.hh \
nmrEltName.hh \
nmrExceptions.hh \
nmrUnit.hh \
partialTokenList.hh \
permutation.hh \
permutationName.hh \
readableNameAssembler.hh 