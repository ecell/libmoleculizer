SUBDIRS = parsers exceptions nodes io

INCLUDES = -I$(top_srcdir) -I../libxml2/include/libxml/ -I../libxml2/include/

DEFS=@DEFS@ -DLIBXMLPP_COMPILATION

h_sources_public = libxml++.h attribute.h dtd.h document.h noncopyable.h keepblanks.h api_export.h
h_sources_private = 

cc_sources = attribute.cc dtd.cc document.cc noncopyable.cc keepblanks.cc
cc_sources_private =

AM_CXXFLAGS=-I../libxml2/include/libxml/

LIBXML2LIB = ../libxml2/libxml2.a

noinst_LTLIBRARIES = libxmlpp.la

# libxmlpp_a_LDFLAGS = @WIN32_LDFLAGS@ -version-info $(GENERIC_LIBRARY_VERSION)
libxmlpp_la_LIBADD = parsers/libparsers.la \
                         exceptions/libexceptions.la \
                         nodes/libnodes.la \
			 io/libio.la \
			../libxml2/libxml2.la

libxmlpp_la_SOURCES = $(cc_sources) $(h_sources_public) $(h_sources_private) $(cc_sources_private)

