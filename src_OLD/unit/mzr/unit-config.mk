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
UNIT_NAME := mzr

# The source files to be compiled, including the .cc extension.  These should
# all be in the cc directory.
CC_FILES := clockDumpable.cc \
	continuator.cc \
	createEvent.cc \
	dumpUtils.cc \
	dumpableNotSpeciesStreamXcpt.cc \
	dumpStateEvent.cc \
	dupDumpableNameXcpt.cc \
	dupSpeciesNameXcpt.cc \
	eventQueue.cc \
	growEvent.cc \
	moleculizer.cc \
	mzrEltName.cc \
	mzrReaction.cc \
	mzrUnit.cc \
	mzrUnitInsert.cc \
	mzrUnitParse.cc \
	noReactEvent.cc \
	noStopEventWarning.cc \
	parametrizer.cc \
	reactCountDumpable.cc \
	reactEventCountDumpable.cc \
	secondsDumpable.cc \
	simTimeDumpable.cc \
	species.cc \
	speciesCountDumpable.cc \
	stopEvent.cc \
	stopEventInPastXcpt.cc \
	tabDumpEvent.cc \
	unitsMgr.cc \
	unkDumpableXcpt.cc \
	unkSpeciesXcpt.cc \
	unkStatStreamXcpt.cc \
	volumeDumpable.cc \
	volumeEvent.cc

# The header files, including the .hh extension.  These should all be in the hh
# directory.
HH_FILES := clockDumpable.hh \
	continuator.hh \
	createEvent.hh \
	dumpableNotSpeciesStreamXcpt.hh \
	dumpStateEvent.hh \
	dumpUtils.hh \
	dupDumpableNameXcpt.hh \
	dupSpeciesNameXcpt.hh \
	eventQueue.hh \
	growEvent.hh \
	inputCapTest.hh \
	inputCapXcpt.hh \
	molarFactor.hh \
	moleculizer.hh \
	mzrEltName.hh \
	mzrEvent.hh \
	mzrReaction.hh \
	mzrSpeciesDumpable.hh \
	mzrSpeciesDumpableImpl.hh \
	mzrSpecies.hh \
	mzrStream.hh \
	mzrUnit.hh \
	noReactEvent.hh \
	noStopEventWarning.hh \
	parametrizer.hh \
	reactCountDumpable.hh \
	reactEventCountDumpable.hh \
	respondReaction.hh \
	secondsDumpable.hh \
	simTimeDumpable.hh \
	speciesCountDumpable.hh \
	stopEvent.hh \
	stopEventInPastXcpt.hh \
	tabDumpEvent.hh \
	unit.hh \
	unitsMgr.hh \
	unkDumpableXcpt.hh \
	unkSpeciesXcpt.hh \
	unkStatStreamXcpt.hh \
	volumeDumpable.hh \
	volumeEvent.hh

# Other units (really libs made by this build) that this
# module requires for linking.
REQUIRED_UNITS := utl

# Libraries that are required for linking this module,
# but are not made by this build.
EXTRA_LIBS :=

# Archive libraries, not made by this build, used in the linking of this
# module.  These locations are typically configured in build-config.mk.
EXTRA_ARCHIVES :=

# Compiler flags and definitions to be used in compiling c++ files for this
# module.
COMPILE_DEFS :=

include $(BUILD)/unit.mk
