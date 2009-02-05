## Process this file with automake to produce Makefile.in.

abs_top_builddir=@abs_top_builddir@

SRCDIR = ..

# One of these needs to be fixed (whatever unit we're writing the unit test for).
LIBFND = ../../fnd/libmoleculizer_fnd.la
LIBUTL = ../../utl/libmoleculizer_utl.la
LIBMZR = ../libmoleculizer_mzr.la
LIBSTOCH = ../../stoch/libmoleculizer_stoch.la
LIBNMR = ../../nmr/libmoleculizer_nmr.la
LIBCPX = ../../cpx/libmoleculizer_cpx.la
LIBMOL = ../../mol/libmoleculizer_mol.la
LIBPLEX = ../../plex/libmoleculizer_plex.la
LIBDIMER = ../../dimer/libmoleculizer_dimer.la
LIBFTR = ../../ftr/libmoleculizer_ftr.la
ALL_MZR_LIBS = $(LIBFND) $(LIBUTL) $(LIBMZR) $(LIBSTOCH) $(LIBNMR) $(LIBNMR) $(LIBCPX) $(LIBMOL) $(LIBPLEX) $(LIBDIMER) $(LIBFTR) 

LIBXML_CFLAGS = @LIBXML_CFLAGS@
LIBXML_LIBS = @LIBXML_LIBS@

GSL_CFLAGS = @GSL_CFLAGS@
GSL_LIBS = @GSL_LIBS@

check_PROGRAMS = 

check_PROGRAMS += moleculizer_test
moleculizer_test_SOURCES = moleculizer_tests.cpp
moleculizer_test_CXXFLAGS= @CXXFLAGS@ -Wall $(LIBXML_CFLAGS)
moleculizer_test_LDADD = $(ALL_MZR_LIBS)

# check_PROGRAMS += c_interface_test
# c_interface_test_SOURCES = c_interface_tests.cpp
# c_interface_test_CXXFLAGS = @CXXFLAGS@ -Wall $(LIBXML_CFLAGS)
# c_interface_test_LDADD = $(ALL_MZR_LIBS)

TESTS=$(check_PROGRAMS)

