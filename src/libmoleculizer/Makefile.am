## Process this file with automake to produce Makefile.in.

SUBDIRS = nauty \
	utl \
	fnd \
	mzr \
	stoch \
	cpx \
	mol \
	nmr \
	plex \
	dimer \
	ftr 

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

lib_LTLIBRARIES = libmoleculizer-1.0.la
libmoleculizer_1_0_la_SOURCES = dummy.cpp
libmoleculizer_1_0_la_LIBADD = \
	$(LIBFND) \
	$(LIBUTL) \
	$(LIBMZR) \
	$(LIBSTOCH) \
	$(LIBNMR) \
	$(LIBCPX) \
	$(LIBMOL) \
	$(LIBPLEX) \
	$(LIBDIMER) \
	$(LIBFTR) \
	$(LIBNAUTY)

