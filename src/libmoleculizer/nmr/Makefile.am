## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_nmr.la

VERSION = @INSTALLATION_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_nmrdir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/nmr

AM_CXXFLAGS = -Wall

libmoleculizer_nmr_la_CXXFLAGS = $(AM_CXXFLAGS) @LIBXMLPP_CFLAGS@ @PYTHON_CSPEC@
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