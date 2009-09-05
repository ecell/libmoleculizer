## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_nmr.la

VERSION = @INSTALLATION_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_nmrdir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/nmr

AM_CXXFLAGS = -Wall

libmoleculizer_nmr_la_CXXFLAGS = $(AM_CXXFLAGS) @LIBXMLPP_CFLAGS@ @PYTHON_CSPEC@
libmoleculizer_nmr_la_SOURCES =\
basicNameAssembler.cpp \
complexOutputState.cpp \
complexSpecies.cpp \
complexSpeciesEncoderNames.cpp \
complexSpeciesOutputMinimizer.cpp \
mangledNameAssembler.cpp \
namedMolecule.cpp \
nmrEltName.cpp \
nmrUnit.cpp \
partialTokenList.cpp \
permutation.cpp \
permutationName.cpp \
readableNameAssembler.cpp 

libmoleculizer_nmr_HEADERS=\
abstractMol.hpp \
basicNameAssembler.hpp \
complexOutputState.hpp \
complexSpeciesEncoderNames.hpp \
complexSpecies.hpp \
complexSpeciesOutputMinimizer.hpp \
mangledNameAssembler.hpp \
nameAssembler.hpp \
nameAssemblers.hpp \
namedMolecule.hpp \
nameEncoderFactory.hpp \
nmrEltName.hpp \
nmrExceptions.hpp \
nmrUnit.hpp \
partialTokenList.hpp \
permutation.hpp \
permutationName.hpp \
readableNameAssembler.hpp 