## Process this file with automake to produce Makefile.in.

SUBDIRS = doc 

if BUILD_LIBXMLPP
SUBDIRS += libxml++-1.0.5
endif

SUBDIRS += src

if ENABLE_DEMOS
SUBDIRS += demos
endif

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libmoleculizer-1.0.pc


