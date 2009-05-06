abs_top_srcdir=@abs_top_srcdir@
SRC=$(abs_top_srcdir)/src
DEMO_SIMULATION_DIR=.

STOCHSIMDIR=$(DEMO_SIMULATION_DIR)/stochastic/src
PARTICLESIMDIR=$(DEMO_SIMULATION_DIR)/particle/src
NETEXPANDDIR=$(DEMO_SIMULATION_DIR)/network_expander/src
CINTERFACEDIR=$(DEMO_SIMULATION_DIR)/c_interface/src
RULEPARSERDIR=$(DEMO_SIMULATION_DIR)/rule_parser/src

bin_PROGRAMS=network_expander_demo stochasticsim_demo particlesim_demo c_interface_demo rule_parser_demo

LIBMZR = $(SRC)/libmoleculizer-1.0.la

LIBXMLPP_CFLAGS = @LIBXMLPP_CFLAGS@
LIBXMLPP_LIBS = @LIBXMLPP_LIBS@

AM_CXXFLAGS = $(LIBXMLPP_CFLAGS) $(GSL_CFLAGS) -Icommon @PYTHON_CPPFLAGS@ -I$(SRC)

noinst_LIBRARIES = libsimbase.a 

libsimbase_a_SOURCES=common/demosimulator.cpp common/demosimulator.hpp
# libsimbase_a_HEADERS=

stochasticsim_demo_SOURCES=\
	$(STOCHSIMDIR)/stochasticsim_main.cpp \
	$(STOCHSIMDIR)/demostochasticsimulator.cpp \
	$(STOCHSIMDIR)/demostochasticsimulator.hpp	

# stochasticsim_demo_HEADERS=

stochasticsim_demo_LDADD = $(LIBMZR) libsimbase.a $(LIBXMLPP_LIBS)

particlesim_demo_SOURCES=\
	$(PARTICLESIMDIR)/particlesim_main.cpp \
	$(PARTICLESIMDIR)/demoparticlesimulator.cpp \
	$(PARTICLESIMDIR)/demoparticlesimulator.hpp

# particlesim_demo_HEADERS=

particlesim_demo_LDADD = $(LIBMZR) libsimbase.a $(LIBXMLPP_LIBS)

network_expander_demo_SOURCES=$(NETEXPANDDIR)/main.cpp
network_expander_demo_LDADD = $(LIBMZR) $(LIBXMLPP_LIBS)

c_interface_demo_SOURCES = $(CINTERFACEDIR)/c_interface_main.c
c_interface_demo_LDADD = $(LIBMZR) $(LIBXMLPP_LIBS)


rule_parser_demo_SOURCES = $(RULEPARSERDIR)/ruleparsermain.cpp
rule_parser_demo_LDADD = $(LIBMZR) $(LIBXMLPP_LIBS)
