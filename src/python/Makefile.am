lib_LTLIBRARIES = libmoleculizer_python.la


TOPSRCDIR=${top_srcdir}
SRCDIR = ..

VERSION = @PACKAGE_VERSION@
MOLECULIZER_DIRNAME = moleculizer-$(VERSION)
libmoleculizer_pythondir=$(includedir)/$(MOLECULIZER_DIRNAME)/python

LIBFND = ../fnd/libmoleculizer_fnd.la
LIBUTL = ../utl/libmoleculizer_utl.la
LIBMZR = ../mzr/libmoleculizer_mzr.la
LIBSTOCH = ../stoch/libmoleculizer_stoch.la
LIBNMR = ../nmr/libmoleculizer_nmr.la
LIBCPX = ../cpx/libmoleculizer_cpx.la
LIBMOL = ../mol/libmoleculizer_mol.la
LIBPLEX = ../plex/libmoleculizer_plex.la
LIBDIMER = ../dimer/libmoleculizer_dimer.la
LIBFTR = ../ftr/libmoleculizer_ftr.la
LIBALL = ../libmoleculizer_super.la
LIBMZRPYTHON = ../python/libmoleculizer_python.la
ALL_MZR_LIBS = $(LIBFND) $(LIBUTL) $(LIBMZR) $(LIBSTOCH) $(LIBNMR) $(LIBNMR) $(LIBCPX) $(LIBMOL) $(LIBPLEX) $(LIBDIMER) $(LIBFTR)

LIBXML_CFLAGS = @LIBXML_CFLAGS@
LIBXML_LIBS = @LIBXML_LIBS@

GSL_CFLAGS = @GSL_CFLAGS@
GSL_LIBS = @GSL_LIBS@

PYTHON=@PYTHON@
PYTHON_VERSION=@PYTHON_VERSION@

LIBPYTHON=-lpython$(PYTHON_VERSION)
LIBBOOSTPYTHON=-lboost_python 

AM_CXXFLAGS=@CXXFLAGS@ -Wall $(LIBXML_CFLAGS) $(GSL_CFLAGS) $(PYTHON_INCLUDES) 
AM_LDFLAGS=$(LIBBOOSTPYTHON)

libmoleculizer_python_la_LIBADD = \
	$(LIBFND) \
	$(LIBUTL) \
	$(LIBMZR) \
	$(LIBSTOCH) \
	$(LIBNMR) \
	$(LIBNMR) \
	$(LIBCPX) \
	$(LIBMOL) \
	$(LIBPLEX) \
	$(LIBDIMER) \
	$(LIBFTR)

libmoleculizer_python_la_SOURCES = 
	MoleculizerPythonWrapper.cpp \
	ReactionNetworkGenerator.cpp



python:
	./make-from-static.sh
