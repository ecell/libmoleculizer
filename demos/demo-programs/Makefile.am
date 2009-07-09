abs_top_srcdir=@abs_top_srcdir@
SRC=$(abs_top_srcdir)/src
DEMO_SIMULATION_DIR=.

bin_PROGRAMS=\
network_expander_demo \
stochasticsim_demo \
particlesim_demo \
c_interface_demo 

LIBMZR = ../../src/libmoleculizer/libmoleculizer-1.0.la

LIBXMLPP_CFLAGS = 
LIBXMLPP_LIBS = @LIBXMLPP_LIBS@

AM_CXXFLAGS = @LIBXMLPP_CFLAGS@ @PYTHON_CSPEC@ -I$(SRC) -Icommon

noinst_LTLIBRARIES = libsimbase.la
libsimbase_la_SOURCES=\
	common/demosimulator.cpp \
	common/demosimulator.hpp

stochasticsim_demo_SOURCES=\
	stochastic/stochasticsim_main.cpp \
	stochastic/demostochasticsimulator.cpp \
	stochastic/demostochasticsimulator.hpp	

stochasticsim_demo_LDADD = libsimbase.la $(LIBMZR) $(LIBXMLPP_LIBS)

particlesim_demo_SOURCES=\
	particle/particlesim_main.cpp \
	particle/demoparticlesimulator.cpp \
	particle/demoparticlesimulator.hpp

particlesim_demo_LDADD = libsimbase.la $(LIBMZR) $(LIBXMLPP_LIBS)

network_expander_demo_SOURCES=network_expander/main.cpp
network_expander_demo_LDADD = $(LIBMZR) $(LIBXMLPP_LIBS)

c_interface_demo_SOURCES = c_interface/c_interface_main.c
c_interface_demo_LDADD = $(LIBMZR) $(LIBXMLPP_LIBS)


