VERSION = @INSTALLATION_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_C_ladir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer


lib_LTLIBRARIES=libmoleculizer_C.la

AM_CXXFLAGS="-I../../libmoleculizer -I../.." @LIBXMLPP_CFLAGS@ @PYTHON_CSPEC@
AM_LDFLAGS = @PYTHON_LSPEC@

libmoleculizer_C_la_SOURCES=\
libmzr_c_interface.cpp	

libmoleculizer_C_la_HEADERS=\
libmzr_c_interface.h