###############################################################################
# Moleculizer - a stochastic simulator for cellular chemistry.
# Copyright (C) 2001  Walter Lawrence (Larry) Lok.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#    
# Contact information:
#   Larry Lok, Research Fellow          Voice: 510-981-8740
#   The Molecular Sciences Institute      Fax: 510-647-0699
#   2168 Shattuck Ave.                  Email: lok@molsci.org
#   Berkeley, CA 94704
###############################################################################

# The name of this module.
UNIT_NAME := cpt

# The source files to be compiled, including the .cc extension.  These should
# all be in the cc directory.
CC_FILES := boundary.cc \
	clockDumpable.cc \
	compartment.cc \
	compartmentGraph.cc \
	compartmentSpecies.cc \
	cptApp.cc \
	cptEltName.cc \
	cptReaction.cc \
	cptUnit.cc \
	cptUnitInsert.cc \
	cptUnitParse.cc \
	createEvent.cc \
	dumpStateEvent.cc \
	dumpableNotSpeciesStreamXcpt.cc \
	dupCompartmentNameXcpt.cc \
	dupConnectionXcpt.cc \
	dupDumpableNameXcpt.cc \
	dupSpeciesNameXcpt.cc \
	globalReaction.cc \
	globalSpecies.cc \
	growEvent.cc \
	noReactEvent.cc \
	noStopEventWarning.cc \
	parseCompartmentGraph.cc \
	parseCreateEvent.cc \
	parseDumpStateEvent.cc \
	parseDumpStream.cc \
	parseGlobalReaction.cc \
	parseGlobalSpecies.cc \
	parseStopEvent.cc \
	propensityDistro.cc \
	reactCountDumpable.cc \
	reactEventCountDumpable.cc \
	secondsDumpable.cc \
	simTimeDumpable.cc \
	singleGlobalSpeciesDumpable.cc \
	speciesCountDumpable.cc \
	stopEvent.cc \
	stopEventInPastXcpt.cc \
	strayEventContentXcpt.cc \
	strayModelContentXcpt.cc \
	strayExplicitSpeciesContentXcpt.cc \
	strayReactionGenXcpt.cc \
	straySpeciesStreamContentXcpt.cc \
	tabDumpEvent.cc \
	unitsMgr.cc \
	unkBoundaryXcpt.cc \
	unkCompartmentXcpt.cc \
	unkDumpableXcpt.cc \
	unkGlobalDumpableXcpt.cc \
	unkSpeciesXcpt.cc \
	unkStatStreamXcpt.cc \
	volumeDumpable.cc \
	volumeEvent.cc

# The header files, including the .hh extension.  These should all be in the hh
# directory.
HH_FILES := boundary.hh \
	clockDumpable.hh \
	compartmentGraph.hh \
	compartment.hh \
	compartmentSpecies.hh \
	cptApp.hh \
	cptEltName.hh \
	cptEvent.hh \
	cptReaction.hh \
	cptUnit.hh \
	createEvent.hh \
	dumpableNotSpeciesStreamXcpt.hh \
	dumpStateEvent.hh \
	dupCompartmentNameXcpt.hh \
	dupConnectionXcpt.hh \
	dupDumpableNameXcpt.hh \
	dupSpeciesNameXcpt.hh \
	eventQueue.hh \
	globalDumpArg.hh \
	globalReaction.hh \
	globalSpeciesDumpableAux.hh \
	globalSpecies.hh \
	growEvent.hh \
	inputCapabilities.hh \
	inputCapTest.hh \
	multiGlobalSpeciesDumpable.hh \
	multiGlobalSpeciesDumpableImpl.hh \
	noReactEvent.hh \
	noStopEventWarning.hh \
	parseCompartmentGraph.hh \
	parseCreateEvent.hh \
	parseDumpStateEvent.hh \
	parseDumpStream.hh \
	parseGlobalReaction.hh \
	parseGlobalSpecies.hh \
	parseStopEvent.hh \
	propensityDistro.hh \
	queryGlobalSpeciesDumpable.hh \
	reactCountDumpable.hh \
	reactEventCountDumpable.hh \
	respondReaction.hh \
	secondsDumpable.hh \
	simTimeDumpable.hh \
	singleGlobalSpeciesDumpable.hh \
	speciesCountDumpable.hh \
	stopEvent.hh \
	stopEventInPastXcpt.hh \
	strayEventContentXcpt.hh \
	strayExplicitSpeciesContentXcpt.hh \
	strayModelContentXcpt.hh \
	strayReactionGenXcpt.hh \
	straySpeciesStreamContentXcpt.hh \
	stream.hh \
	tabDumpEvent.hh \
	unit.hh \
	unitsMgr.hh \
	unkBoundaryXcpt.hh \
	unkCompartmentXcpt.hh \
	unkDumpableXcpt.hh \
	unkGlobalDumpableXcpt.hh \
	unkSpeciesXcpt.hh \
	unkStatStreamXcpt.hh \
	volumeDumpable.hh \
	volumeEvent.hh

# Other units (really libs made by this build) that this
# module requires for linking.
REQUIRED_UNITS := utl

# Libraries that are required for linking this module,
# but are not made by this build.
EXTRA_LIBS := $(GSL_LIB) \
	$(GSL_CBLAS_LIB)

# Archive libraries, not made by this build, used in the linking of this
# module.  These locations are typically configured in build-config.mk.
EXTRA_ARCHIVES :=

# Compiler flags and definitions to be used in compiling c++ files for this
# module.
COMPILE_DEFS :=

include $(BUILD)/unit.mk
