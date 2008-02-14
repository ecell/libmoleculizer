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

# This makefile handles the simplest transformations of schema files by xslt:
# 1. Removal of element/attribute documentation elements.
# 2. Generation of the compact form of the rng schemas.
# Flattening is more complex and I plan to give it its own makefile.

# Where schemas containing documentation come from.
DOC_SCHEMA_SOURCE_DIR := $(SCHEMA)/schema-doc

PREEN_LIST += $(DOC_SCHEMA_SOURCE_DIR)/*.bak \
	$(DOC_SCHEMA_SOURCE_DIR)/*~

# Schema-doc schemas.  These need documentation removal if they are to be
# used as standard rng schemas.
DOC_SCHEMAS := bndKinase.rng \
	cdm.rng \
	cft.rng \
	clx.rng \
	cml.rng \
	compartment.rng \
	cst.rng \
	dimer.rng \
	ftr.rng \
	gpa.rng \
	modKinase.rng \
	moleculizer.rng \
	mol.rng \
	mzrState.rng \
	nucEx.rng \
	odie.rng \
	plex.rng \
	rk4ode.rng \
	rk4tau.rng \
	scaffold.rng \
	stoch.rng

# Where the fixed schemas that don't need documentation removal come from.
FIXED_SCHEMA_SOURCE_DIR := $(SCHEMA)/fixed

PREEN_LIST += $(FIXED_SCHEMA_SOURCE_DIR)/*.bak \
	$(FIXED_SCHEMA_SOURCE_DIR)/*~

# Fixed schemas that don't need documentation removal.
FIXED_SCHEMAS := mzr-defaults.rng \
	rngdoc.rng

#######################

# Where to put processed schemas, etc. in the build area.
SCHEMA_TARGET_DIR := $(XML)/schema

$(SCHEMA_TARGET_DIR) : | $(XML)
	mkdir $@

# Convenience target to build all schemas.
$(SCHEMA_TARGET_DIR)/.schemas : | $(SCHEMA_TARGET_DIR)
	touch $@

# Add all the schemas to the general xml target.
xml-target : $(SCHEMA_TARGET_DIR)/.schemas

#######################

# Where to put schemas with documentation removed.
DOC_REMOVED_DIR := $(SCHEMA_TARGET_DIR)/doc-removed

$(DOC_REMOVED_DIR) : $(SCHEMA_TARGET_DIR)
	mkdir $@

DOC_REMOVED_SCHEMAS := $(addprefix $(DOC_REMOVED_DIR)/,$(DOC_SCHEMAS))

$(DOC_REMOVED_SCHEMAS) : | $(DOC_REMOVED_DIR)

$(DOC_REMOVED_SCHEMAS) : $(DOC_REMOVED_DIR)/% : $(DOC_SCHEMA_SOURCE_DIR)/%
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/doc2rng.xsl \
	-xml \
	-out $@

# Convenience target to make all schema with doc removed.
$(DOC_REMOVED_DIR)/.doc-removed : $(DOC_REMOVED_SCHEMAS)
	touch $@

# Add the schemas with documentation removed to the general schema target.
$(SCHEMA_TARGET_DIR)/.schemas : $(DOC_REMOVED_DIR)/.doc-removed

#######################

# Compact schemas are created from plain rng schemas, which in turn come from
# removing element/attribute documentation from doc-schemas.

# Where to put compact schemas.
COMPACT_DIR := $(SCHEMA_TARGET_DIR)/compact

$(COMPACT_DIR) : $(SCHEMA_TARGET_DIR)
	mkdir $@

# We also change the file extension to .rnc for compact schemas.
COMPACT_SCHEMAS := $(addprefix $(COMPACT_DIR)/,$(DOC_SCHEMAS:.rng=.rnc))

$(COMPACT_SCHEMAS) : | $(COMPACT_DIR)

$(COMPACT_SCHEMAS) : $(COMPACT_DIR)/%.rnc : $(DOC_REMOVED_DIR)/%.rng
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/RngToRnc-1_4/RngToRncText.xsl \
	-xml \
	-out $@

# Convenience target to make all compact schema.
$(COMPACT_DIR)/.compact : $(COMPACT_SCHEMAS)
	touch $@

# Add the compact schemas to the general schema target.
$(SCHEMA_TARGET_DIR)/.schemas : $(COMPACT_DIR)/.compact

#######################

include $(BUILD)/schema-flat.mk
