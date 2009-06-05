## Process this file with automake to produce Makefile.in.

abs_top_builddir=@abs_top_builddir@


MOLECULIZER_DIRNAME = @MOLECULIZER_DIRNAME@

moleculizerheaderdir=${include_dir}/${MOLECULIZER_DIRNAME}

SRCDIR = $(abs_top_builddir)/src

LIBFND = libmoleculizer/fnd/libmoleculizer_fnd.la
LIBUTL = libmoleculizer/utl/libmoleculizer_utl.la
LIBMZR = libmoleculizer/mzr/libmoleculizer_mzr.la
LIBSTOCH = libmoleculizer/stoch/libmoleculizer_stoch.la
LIBNMR = libmoleculizer/nmr/libmoleculizer_nmr.la
LIBCPX = libmoleculizer/cpx/libmoleculizer_cpx.la
LIBMOL = libmoleculizer/mol/libmoleculizer_mol.la
LIBPLEX = libmoleculizer/plex/libmoleculizer_plex.la
LIBDIMER = libmoleculizer/dimer/libmoleculizer_dimer.la
LIBFTR = libmoleculizer/ftr/libmoleculizer_ftr.la
LIBNAUTY = libmoleculizer/nauty/libmoleculizer_nauty.la

ALL_MZR_LIBS = $(LIBFND) $(LIBUTL) $(LIBMZR) $(LIBSTOCH) $(LIBNMR) $(LIBNMR) $(LIBCPX) $(LIBMOL) $(LIBPLEX) $(LIBDIMER) $(LIBFTR) $(LIBNAUTY)	

moleculizerdir = $(SRCDIR)

LIBXMLPP_CFLAGS = @LIBXMLPP_CFLAGS@
LIBXMLPP_LIBS = @LIBXMLPP_LIBS@

AM_CXXFLAGS = $(LIBXMLPP_CFLAGS)
AM_LDFLAGS = @PYTHON_LDFLAGS@ @PYTHON_EXTRA_LDFLAGS@ 

lib_LTLIBRARIES = libmoleculizer-1.0.la
libmoleculizer_1_0_la_SOURCES =
# nodist_EXTRA_libmoleculizer_1_0_la_SOURCES = dummy.cpp

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


SUBDIRS = libmoleculizer/nauty \
	libmoleculizer/utl \
	libmoleculizer/fnd \
	libmoleculizer/mzr \
	libmoleculizer/stoch \
	libmoleculizer/cpx \
	libmoleculizer/mol \
	libmoleculizer/nmr \
	libmoleculizer/plex \
	libmoleculizer/dimer \
	libmoleculizer/ftr 

# SUBDIRS += python