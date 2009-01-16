abs_top_builddir=@abs_top_builddir@

SRC=$(abs_top_builddir)/src
DEMO_SIMULATION_DIR=$(abs_top_builddir)/doc/demos/simulators

STOCHSIMDIR=$(DEMO_SIMULATION_DIR)/stochastic/src
PARTICLESIMDIR=$(DEMO_SIMULATION_DIR)/particle/src
NETEXPANDDIR=$(DEMO_SIMULATION_DIR)/network_expander/src
CINTERFACEDIR=$(DEMO_SIMULATION_DIR)/c_interface/src
CPARTICLEDIR=$(DEMO_SIMULATION_DIR)/c_particle/src

bin_PROGRAMS=network_expander_demo stochasticsim_demo particlesim_demo c_interface_demo 

LIBFND = $(SRC)/fnd/libmoleculizer_fnd.la
LIBUTL = $(SRC)/utl/libmoleculizer_utl.la
LIBMZR = $(SRC)/mzr/libmoleculizer_mzr.la
LIBSTOCH = $(SRC)/stoch/libmoleculizer_stoch.la
LIBNMR = $(SRC)/nmr/libmoleculizer_nmr.la
LIBCPX = $(SRC)/cpx/libmoleculizer_cpx.la
LIBMOL = $(SRC)/mol/libmoleculizer_mol.la
LIBPLEX = $(SRC)/plex/libmoleculizer_plex.la
LIBDIMER = $(SRC)/dimer/libmoleculizer_dimer.la
LIBFTR = $(SRC)/ftr/libmoleculizer_ftr.la

ALL_MZR_LIBS = $(LIBFND) $(LIBUTL) $(LIBMZR) $(LIBSTOCH) $(LIBNMR) $(LIBNMR) $(LIBCPX) $(LIBMOL) $(LIBPLEX) $(LIBDIMER) $(LIBFTR) 

LIBXML_CFLAGS = @LIBXML_CFLAGS@
LIBXML_LIBS = @LIBXML_LIBS@

AM_CXXFLAGS = $(LIBXML_CFLAGS) $(GSL_CFLAGS) -Icommon

noinst_LIBRARIES = libsimbase.a 

libsimbase_a_SOURCES=common/demosimulator.cpp

stochasticsim_demo_SOURCES=\
	$(STOCHSIMDIR)/stochasticsim_main.cpp \
	$(STOCHSIMDIR)/demostochasticsimulator.cpp

stochasticsim_demo_LDADD = $(ALL_MZR_LIBS) libsimbase.a

particlesim_demo_SOURCES=\
	$(PARTICLESIMDIR)/particlesim_main.cpp \
	$(PARTICLESIMDIR)/demoparticlesimulator.cpp

particlesim_demo_LDADD = $(ALL_MZR_LIBS) libsimbase.a

network_expander_demo_SOURCES=$(NETEXPANDDIR)/main.cpp
network_expander_demo_LDADD = $(ALL_MZR_LIBS)

c_interface_demo_SOURCES = $(CINTERFACEDIR)/c_interface_main.c
c_interface_demo_LDADD = $(ALL_MZR_LIBS)

# c_particle_demo_SOURCES = $(CPARTICLEDIR)/cparticle_main.c
# c_particle_demo_LDADD = $(ALL_MZR_LIBS)