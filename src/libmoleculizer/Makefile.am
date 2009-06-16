## Process this file with automake to produce Makefile.in.


LIBXMLPP_CFLAGS = @LIBXMLPP_CFLAGS@
PYTHON_CFLAGS = @PYTHON_CPPFLAGS@

LIBXMLPP_LIBS = @LIBXMLPP_LIBS@


AM_CXXFLAGS = $(LIBXMLPP_CFLAGS)

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

lib_LTLIBRARIES = libmoleculizer-1.0.la

libmoleculizer_1_0_la_SOURCES =

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



if BUILD_LIBXMLPP
libmoleculizer_1_0_la_LIBADD += ../../libxml++-1.0.5/libxml++/libxmlpp.a
endif



EXTRA_DIST=Makefile.am