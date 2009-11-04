abs_top_srcdir=@abs_top_srcdir@
SRC=$(abs_top_srcdir)/src
DEMO_SIMULATION_DIR=.

bin_PROGRAMS=\
mzrnetworkgenerator \
ssa_integration_demo \
particle_sim_integration_demo \
c_interface_demo 

LIBMZR = ../../src/libmoleculizer/libmoleculizer-1.0.la
LIBCMZR = ../../src/bindings/C/libmoleculizer_C.la

LIBXMLPP_CFLAGS = 
LIBXMLPP_LIBS = @LIBXMLPP_LIBS@

AM_CXXFLAGS = @LIBXMLPP_CFLAGS@ @PYTHON_CSPEC@ -I$(SRC) -Icommon -I$(SRC)/bindings/C @LIBSBML_CFLAGS@

noinst_LTLIBRARIES = libsimbase.la
libsimbase_la_SOURCES=\
	common/demosimulator.cpp \
	common/demosimulator.hpp

ssa_integration_demo_SOURCES=\
	ssa_integration_demo_src/stochasticsim_main.cpp \
	ssa_integration_demo_src/demostochasticsimulator.cpp \
	ssa_integration_demo_src/demostochasticsimulator.hpp	

ssa_integration_demo_LDADD = libsimbase.la $(LIBMZR) $(LIBXMLPP_LIBS)

particle_sim_integration_demo_SOURCES=\
	particle_sim_integration_demo_src/particlesim_main.cpp \
	particle_sim_integration_demo_src/demoparticlesimulator.cpp \
	particle_sim_integration_demo_src/demoparticlesimulator.hpp

particle_sim_integration_demo_LDADD = libsimbase.la $(LIBMZR) $(LIBXMLPP_LIBS)

mzrnetworkgenerator_SOURCES=mzrnetworkgenerator_src/mzrnetgen.cpp
mzrnetworkgenerator_LDADD = $(LIBMZR) $(LIBXMLPP_LIBS) @LIBSBML_LIBS@

AM_CFLAGS= -I../bindings/C

c_interface_demo_SOURCES = c_interface_demo_src/c_interface_main.c
c_interface_demo_LDADD = $(LIBCMZR) $(LIBMZR) $(LIBXMLPP_LIBS)


