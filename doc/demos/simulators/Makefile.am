

SRC=../../../src
DEMO_SIMULATION_DIR=.

STOCHSIMDIR=$(DEMO_SIMULATION_DIR)/stochastic/src
PARTICLESIMDIR=$(DEMO_SIMULATION_DIR)/particle/src
NETEXPANDDIR=$(DEMO_SIMULATION_DIR)/network_expander/src
CINTERFACEDIR=$(DEMO_SIMULATION_DIR)/c_interface/src

bin_PROGRAMS=network_expander_demo stochasticsim_demo particlesim_demo c_interface_demo 

LIBMZR = $(SRC)/libmoleculizer-1.0.la

LIBXML_CFLAGS = @LIBXML_CFLAGS@
LIBXML_LIBS = @LIBXML_LIBS@

AM_CXXFLAGS = $(LIBXML_CFLAGS) $(GSL_CFLAGS) -Icommon

noinst_LIBRARIES = libsimbase.a 

libsimbase_a_SOURCES=common/demosimulator.cpp

stochasticsim_demo_SOURCES=\
	$(STOCHSIMDIR)/stochasticsim_main.cpp \
	$(STOCHSIMDIR)/demostochasticsimulator.cpp

stochasticsim_demo_LDADD = $(LIBMZR) libsimbase.a $(LIBXML_LIBS)

particlesim_demo_SOURCES=\
	$(PARTICLESIMDIR)/particlesim_main.cpp \
	$(PARTICLESIMDIR)/demoparticlesimulator.cpp

particlesim_demo_LDADD = $(LIBMZR) libsimbase.a $(LIBXML_LIBS)

network_expander_demo_SOURCES=$(NETEXPANDDIR)/main.cpp
network_expander_demo_LDADD = $(LIBMZR) $(LIBXML_LIBS)

c_interface_demo_SOURCES = $(CINTERFACEDIR)/c_interface_main.c
c_interface_demo_LDADD = $(LIBMZR) $(LIBXML_LIBS)
