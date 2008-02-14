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
UNIT_NAME := cpx

# The source files to be compiled, including the .cc extension.  These should
# all be in the cc directory.
CC_FILES := badModMolStateXcpt.cc \
	basicBndSite.cc \
	modMolMixin.cc \
	modPlexQueryTypeXcpt.cc \
	modStateMixin.cc \
	noKineticConstsXcpt.cc \
	plexIso.cc \
	plexMap.cc \
	plexNotConnectedXcpt.cc \
	siteToShapeMap.cc \
	unmappedSiteSpecXcpt.cc

# The header files, including the .hh extension.  These should all be in the hh
# directory.
HH_FILES := alloMol.hh \
	alloMolImpl.hh \
	badModMolStateXcpt.hh \
	basicBndSite.hh \
	basicMol.hh \
	basicMolImpl.hh \
	basicPlex.hh \
	basicPlexImpl.hh \
	binding.hh \
	cxBinding.hh \
	cxBindingImpl.hh \
	cxMol.hh \
	cxMolImpl.hh \
	cxOmni.hh \
	cxOmniImpl.hh \
	cxSite.hh \
	cxSiteImpl.hh \
	ftrSpec.hh \
	hashMolRec.hh \
	hashMolRecImpl.hh \
	isoSearch.hh \
	isoSearchImpl.hh \
	knownBindings.hh \
	modification.hh \
	modMixinQuery.hh \
	modMol.hh \
	modMolImpl.hh \
	modMolMixin.hh \
	modMolState.hh \
	modMolStateQuery.hh \
	modPlexQueryTypeXcpt.hh \
	modStateMixin.hh \
	molState.hh \
	molStateQuery.hh \
	noKineticConstsXcpt.hh \
	omniPlexFeature.hh \
	omniPlex.hh \
	omniStructureQuery.hh \
	omniStructureQueryImpl.hh \
	plexFamily.hh \
	plexFamilyImpl.hh \
	plexIso.hh \
	plexIsoImpl.hh \
	plexMap.hh \
	plexMapImpl.hh \
	plexNotConnectedXcpt.hh \
	plexQuery.hh \
	plexQueryImpl.hh \
	plexSpcsMixin.hh \
	plexSpcsMixinImpl.hh \
	prm.hh \
	queryAlloList.hh \
	queryAlloListImpl.hh \
	recognizer.hh \
	recognizerImpl.hh \
	reportIsoSearch.hh \
	siteShape.hh \
	siteToShapeMap.hh \
	siteToShapeMapImpl.hh \
	smallMol.hh \
	stateMol.hh \
	structuralBinding.hh \
	subPlexSpec.hh \
	unmappedSiteSpecXcpt.hh

# Other units (really libs made by this build) that this
# module requires for linking.
REQUIRED_UNITS :=

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
