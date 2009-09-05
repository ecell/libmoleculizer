## Process this file with automake to produce Makefile.in.

noinst_LTLIBRARIES = libmoleculizer_fnd.la

VERSION = @INSTALLATION_VERSION@
MOLECULIZER_DIRNAME = libmoleculizer-$(VERSION)
libmoleculizer_fnddir=$(includedir)/$(MOLECULIZER_DIRNAME)/libmoleculizer/fnd

AM_CXXFLAGS = -Wall

libmoleculizer_fnd_la_SOURCES =\
dumpStream.cpp \
fndXcpt.cpp \
pchem.cpp \
physConst.cpp

libmoleculizer_fnd_HEADERS=\
basicDmpColumn.hpp \
basicDumpable.hpp \
basicReaction.hpp \
basicSpecies.hpp \
binaryRxnGen.hpp \
coreRxnGen.hpp \
dmpColumn.hpp \
dumpStream.hpp \
dumpable.hpp \
event.hpp \
feature.hpp \
featureContext.hpp \
featureMap.hpp \
featureStimulus.hpp \
fndXcpt.hpp \
gillspReaction.hpp \
massive.hpp \
multiSpeciesDumpable.hpp \
newContextStimulus.hpp \
newSpeciesStimulus.hpp \
notifier.hpp \
pchem.hpp \
physConst.hpp \
query.hpp \
queryImpl.hpp \
querySpeciesDumpable.hpp \
reactionNetworkDescription.hpp \
reactionNetworkDescriptionImpl.hpp \
reactionNetworkComponent.hpp \
rxnGen.hpp \
sensitive.hpp \
sensitiveOnce.hpp \
sensitivityList.hpp \
stateVar.hpp \
varDumpable.hpp

