## Process this file with automake to produce Makefile.in.

lib_LTLIBRARIES = libmoleculizer_nmr.la

abs_top_builddir=@abs_top_builddir@

## SRCDIR = $(abs_top_builddir)/src
SRCDIR = ..

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

if STATIC_BUILD
AM_LDFLAGS= -all-static
endif

libmoleculizer_nmr_la_SOURCES =\
ComplexSpeciesManglerNames.cc \
newMol.cc \
nmrEltName.cc \
permutation.cc

libmoleculizer_nmr_HEADERS=\
abstractMol.hh \
basicNameAssembler.hh \
complexOutputState.hh \
complexSpecies.hh \
complexSpeciesImpl.hh \
complexSpeciesOutputMinimizer.hh \
complexSpeciesOutputMinimizerImpl.hh \
ComplexSpeciesManglerNames.hh \
csException.hh \
csUtl.hh \
mangledNameAssembler.hh \
mangledNameAssemblerImpl.hh \
nameAssembler.hh \
newMol.hh \
nmrEltName.hh \
partialTokenList.hh \
permutation.hh \
permutationName.hh \
readableNameAssembler.hh

