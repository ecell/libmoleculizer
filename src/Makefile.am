# ./src/Makefile.am
abs_top_builddir=@abs_top_builddir@

# This is the order in which the libraries should be built.
# (dependencies go from left to right)
SUBDIRS = fnd utl mzr stoch cpx mol nmr plex dimer ftr

# This is the name ./moleculizer-2.0, for the purpose of having
# /usr/local/include/moleculizer-2.0/ and possibly something 
# like /usr/local/lib/moleculizer-2.0/.

MOLECULIZER_DIRNAME = @MOLECULIZER_DIRNAME@
moleculizerheaderdir=${include_dir}/${MOLECULIZER_DIRNAME}

SRCDIR = $(abs_top_builddir)/src

# This is the local target directory name.
# I think it has something to do with header dependencies.
moleculizerdir = $(SRCDIR)

# These are constants defining the different library
# dependencies (clearly).

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

# This lib add, with this syntax (I think) was taken from ecell/libecs/Makefile.am

bin_PROGRAMS = moleculizer

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
$(LIBFTR)

AM_CXXFLAGS = -I/usr/include/libxml++-1.0

moleculizer_DEPENDENCIES = $(ALL_MZR_LIBS)

moleculizer_SOURCES=main.cc

all-local: moleculizer
