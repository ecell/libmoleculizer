## Process this file with automake to produce Makefile.in.

abs_top_builddir=@abs_top_builddir@

SUBDIRS = utl \
	fnd \
	mzr \
	stoch \
	cpx \
	mol \
	nmr \
	plex \
	dimer \
	ftr 
## 	python 

MOLECULIZER_DIRNAME = @MOLECULIZER_DIRNAME@
moleculizerheaderdir=${include_dir}/${MOLECULIZER_DIRNAME}

SRCDIR = $(abs_top_builddir)/src

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
ALL_MZR_LIBS = $(LIBFND) $(LIBUTL) $(LIBMZR) $(LIBSTOCH) $(LIBNMR) $(LIBNMR) $(LIBCPX) $(LIBMOL) $(LIBPLEX) $(LIBDIMER) $(LIBFTR)

## This is the local target directory name.
## I think it has something to do with header dependencies.
moleculizerdir = $(SRCDIR)

LIBXML_CFLAGS = @LIBXML_CFLAGS@
LIBXML_LIBS = @LIBXML_LIBS@

AM_CXXFLAGS = $(LIBXML_CFLAGS)
AM_CFLAGS=-g -O0

lib_LTLIBRARIES = libmoleculizer.la
libmoleculizer_la_SOURCES =
nodist_EXTRA_libmoleculizer_la_SOURCES = dummy.cpp

libmoleculizer_la_LIBADD = \
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

# moleculizer_LDADD = \
# 	$(LIBFND) \
# 	$(LIBUTL) \
# 	$(LIBMZR) \
# 	$(LIBSTOCH) \
# 	$(LIBNMR) \
# 	$(LIBCPX) \
# 	$(LIBMOL) \
# 	$(LIBPLEX) \
# 	$(LIBDIMER) \
# 	$(LIBFTR)


bin_PROGRAMS = moleculizer 

if STATIC_BUILD
AM_LDFLAGS=-all-static
endif

moleculizer_DEPENDENCIES = libmoleculizer.la

moleculizer_SOURCES=main.cc
moleculizer_LDADD = libmoleculizer.la \
	$(LIBXML_LIBS)



