## Process this file with automake to produce Makefile.in.

SUBDIRS = nauty utl fnd mzr stoch cpx mol nmr plex dimer ftr 

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
LIBNAUTY = nauty/libmoleculizer_nauty.la

lib_LTLIBRARIES = libmoleculizer-1.0.la

libmoleculizer_1_0_la_SOURCES = dummy.cpp
libmoleculizer_1_0_la_LDFLAGS = @PYTHON_LSPEC@
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

