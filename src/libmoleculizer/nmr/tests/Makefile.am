## Process this file with automake to produce Makefile.in.

abs_top_builddir=@abs_top_builddir@

SRCDIR = ..

VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = moleculizer-$(VERSION)

LIBFND = ../../fnd/libmoleculizer_fnd.la
LIBUTL = ../../utl/libmoleculizer_utl.la
LIBMZR = ../../mzr/libmoleculizer_mzr.la
LIBSTOCH = ../../stoch/libmoleculizer_stoch.la
LIBNMR = ../libmoleculizer_nmr.la
LIBCPX = ../../cpx/libmoleculizer_cpx.la
LIBMOL = ../../mol/libmoleculizer_mol.la
LIBPLEX = ../../plex/libmoleculizer_plex.la
LIBDIMER = ../../dimer/libmoleculizer_dimer.la
LIBFTR = ../../ftr/libmoleculizer_ftr.la
ALL_MZR_LIBS = $(LIBFND) $(LIBUTL) $(LIBMZR) $(LIBSTOCH) $(LIBNMR) $(LIBNMR) $(LIBCPX) $(LIBMOL) $(LIBPLEX) $(LIBDIMER) $(LIBFTR)

check_PROGRAMS = 

check_PROGRAMS += permutation_test
permutation_test_SOURCES = permutation_tests.cpp
permutation_test_LDADD = $(ALL_MZR_LIBS)

check_PROGRAMS += mol_test
mol_test_SOURCES = mol_tests.cpp
mol_test_LDADD = $(ALL_MZR_LIBS)

check_PROGRAMS += complex_species_test
complex_species_test_SOURCES = complexSpecies_test.cpp
complex_species_test_LDADD = $(ALL_MZR_LIBS)

# check_PROGRAMS += namer_test
# namer_test_SOURCES =namer_tests.cpp
# namer_test_LDADD = $(ALL_MZR_LIBS)

TESTS=$(check_PROGRAMS)

