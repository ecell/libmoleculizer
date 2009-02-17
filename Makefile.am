## Process this file with automake to produce Makefile.in.

abs_top_builddir=@abs_top_builddir@

XMLOPERATOR_DIR=$(abs_top_builddir)/xmloperator_2_3_0


SUBDIRS = src doc

if ENABLE_DEMOS
SUBDIRS += demos
endif

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libmoleculizer-1.0.pc


