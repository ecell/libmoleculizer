## Process this file with automake to produce Makefile.in.

abs_top_builddir=@abs_top_builddir@

SRCDIR = ..

# One of these needs to be fixed (whatever unit we're writing the unit test for).
LIBFND = ../../fnd/libmoleculizer_fnd.la
LIBUTL = ../../utl/libmoleculizer_utl.la
LIBMZR = ../../mzr/libmoleculizer_mzr.la
LIBSTOCH = ../../stoch/libmoleculizer_stoch.la
LIBNMR = ../../nmr/libmoleculizer_nmr.la
LIBCPX = ../../cpx/libmoleculizer_cpx.la
LIBMOL = ../../mol/libmoleculizer_mol.la
LIBPLEX = ../../plex/libmoleculizer_plex.la
LIBDIMER = ../../dimer/libmoleculizer_dimer.la
LIBFTR = ../../ftr/libmoleculizer_ftr.la
ALL_MZR_LIBS = $(LIBFND) $(LIBUTL) $(LIBMZR) $(LIBSTOCH) $(LIBNMR) $(LIBNMR) $(LIBCPX) $(LIBMOL) $(LIBPLEX) $(LIBDIMER) $(LIBFTR)

check_PROGRAMS = 

check_PROGRAMS += template_test
template_test_SOURCES =
template_test_LDADD = $(ALL_MZR_LIBS)


TESTS=$(check_PROGRAMS)

