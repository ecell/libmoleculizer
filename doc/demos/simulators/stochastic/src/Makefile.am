abs_top_builddir=@abs_top_builddir@
SRC=$(abs_top_builddir)/src

bin_PROGRAMS=stochasticsim_demo

LIBFND = $(SRC)/fnd/libmoleculizer_fnd.la
LIBUTL = $(SRC)/utl/libmoleculizer_utl.la
LIBMZR = $(SRC)/mzr/libmoleculizer_mzr.la
LIBSTOCH = $(SRC)/stoch/libmoleculizer_stoch.la
LIBNMR = $(SRC)/nmr/libmoleculizer_nmr.la
LIBCPX = $(SRC)/cpx/libmoleculizer_cpx.la
LIBMOL = $(SRC)/mol/libmoleculizer_mol.la
LIBPLEX = $(SRC)/plex/libmoleculizer_plex.la
LIBDIMER = $(SRC)/dimer/libmoleculizer_dimer.la
LIBFTR = $(SRC)/ftr/libmoleculizer_ftr.la
LIBNAUTY = $(SRC)/nauty/libmoleculizer_nauty.a

ALL_MZR_LIBS = $(LIBFND) $(LIBUTL) $(LIBMZR) $(LIBSTOCH) $(LIBNMR) $(LIBNMR) $(LIBCPX) $(LIBMOL) $(LIBPLEX) $(LIBDIMER) $(LIBFTR) $(LIBNAUTY)

LIBXML_CFLAGS = @LIBXML_CFLAGS@
LIBXML_LIBS = @LIBXML_LIBS@

GSL_CFLAGS = @GSL_CFLAGS@
GSL_LIBS = @GSL_LIBS@

AM_CXXFLAGS = $(LIBXML_CFLAGS) $(GSL_CFLAGS)

stochasticsim_demo_SOURCES=main.cpp \
	demosimulator.cpp \
	demostochasticsimulator.cpp

stochasticsim_demo_LDADD = $(ALL_MZR_LIBS)

