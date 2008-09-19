## Process this file with automake to produce Makefile.in.

abs_top_builddir=@abs_top_builddir@

## This is the order in which the libraries should be built.
## (dependencies go from left to right)
SUBDIRS = doc \
	nauty \
	fnd \
	utl \
	mzr \
	stoch \
	cpx \
	mol \
	nmr \
	plex \
	dimer \
	ftr \
	python 

MOLECULIZER_DIRNAME = @MOLECULIZER_DIRNAME@
moleculizerheaderdir=${include_dir}/${MOLECULIZER_DIRNAME}

SRCDIR = $(abs_top_builddir)/src

## This is the local target directory name.
## I think it has something to do with header dependencies.
moleculizerdir = $(SRCDIR)

## These are constants defining the different library
## dependencies (clearly).
LIBFND = fnd/libmoleculizer_fnd.la
LIBUTL = utl/libmoleculizer_utl.la
LIBMZR = mzr/libmoleculizer_mzr.la
LIBSTOCH = stoch/libmoleculizer_stoch.la
LIBNMR = nmr/libmoleculizer_nmr.la
LIBCPX = cpx/libmoleculizer_cpx.la
LIBMOL = mol/libmoleculizer_mol.la
LIBPLEX = plex/libmoleculizer_plex.la
LIBDIMER = dimer/libmoleculizer_dimer.la
LIBFTR = ftr/libmoleculizer_ftr.la
LIBNAUTY = nauty/libmoleculizer_nauty.a

ALL_MZR_LIBS = $(LIBFND) $(LIBUTL) $(LIBMZR) $(LIBSTOCH) $(LIBNMR) $(LIBNMR) $(LIBCPX) $(LIBMOL) $(LIBPLEX) $(LIBDIMER) $(LIBFTR) $(LIBNAUTY)

bin_PROGRAMS = moleculizer

lib_LTLIBRARIES = libmoleculizer_super.la

libmoleculizer_super_la_SOURCES =

nodist_EXTRA_libmoleculizer_super_la_SOURCES = dummy.cpp
libmoleculizer_super_la_LIBADD = \
	$(LIBFND) \
	$(LIBUTL) \
	$(LIBMZR) \
	$(LIBSTOCH) \
	$(LIBNMR) \
	$(LIBCPX) \
	$(LIBMOL) \
	$(LIBPLEX) \
	$(LIBDIMER) \
	$(LIBFTR)

moleculizer_LDADD = \
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

LIBXML_CFLAGS = @LIBXML_CFLAGS@
LIBXML_LIBS = @LIBXML_LIBS@

AM_CXXFLAGS = $(LIBXML_CFLAGS)

if STATIC_BUILD
AM_LDFLAGS=-static
endif

moleculizer_DEPENDENCIES = $(ALL_MZR_LIBS)

moleculizer_SOURCES=main.cc

