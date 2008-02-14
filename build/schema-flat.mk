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

# This makefile shows how to "flatten" doc-schema, in the sense of satisfying
# includes.  This is to make xslt processing of the doc-schema file easier.

# Schema files we want in flattened form.
FLAT_SCHEMA_NAMES := compartment.rng \
	moleculizer.rng \
	odie.rng \
	rk4tau.rng \
	rk4ode.rng \
	mzrState.rng

# Schema files included in moleculizer.rng and mzrState.rng.
MOLECULIZER_SCHEMATA := dimer.rng \
	gpa.rng \
	modKinase.rng \
	mol.rng \
	moleculizer.rng \
	nucEx.rng \
	plex.rng \
	scaffold.rng \
	stoch.rng \
	bndKinase.rng \
	ftr.rng

MZR_SCMTA := $(addprefix $(DOC_SCHEMA_SOURCE_DIR)/,$(MOLECULIZER_SCHEMATA))

# Schema files included in compartment.rng.
CPT_SCHEMATA := cst.rng \
	cml.rng \
	clx.rng \
	cdm.rng \
	cft.rng

CPT_SCMTA := $(addprefix $(DOC_SCHEMA_SOURCE_DIR)/,$(CPT_SCHEMATA))

# Where to put flattened schemas.
FLATTENED_DIR := $(SCHEMA_TARGET_DIR)/flat

$(FLATTENED_DIR) : | $(SCHEMA_TARGET_DIR)
	mkdir $@

################

# A convenience target to make all the flat schemas.

FLATTENED_SCHEMA := $(addprefix $(FLATTENED_DIR)/,$(FLAT_SCHEMA_NAMES))

$(FLATTENED_DIR)/.flattened : $(FLATTENED_SCHEMA)
	touch $@

# Add the flat schemas to the general schema target.
$(SCHEMA_TARGET_DIR)/.schemas : $(FLATTENED_DIR)/.flattened

################

# Note that moleculizer, mzrState, and compartment are slightly different
# from one another and from the rest.

$(FLATTENED_DIR)/moleculizer.rng : | $(FLATTENED_DIR)

$(FLATTENED_DIR)/moleculizer.rng : $(DOC_SCHEMA_SOURCE_DIR)/moleculizer.rng $(MZR_SCMTA)
	java org.apache.xalan.xslt.Process \
	-in $(DOC_SCHEMA_SOURCE_DIR)/moleculizer.rng \
	-xsl $(XSL)/flatten-doc.xsl \
	-xml \
	-out $@

$(FLATTENED_DIR)/mzrState.rng : | $(FLATTENED_DIR)

$(FLATTENED_DIR)/mzrState.rng : $(DOC_SCHEMA_SOURCE_DIR)/mzrState.rng $(MZR_SCMTA)
	java org.apache.xalan.xslt.Process \
	-in $(DOC_SCHEMA_SOURCE_DIR)/mzrState.rng \
	-xsl $(XSL)/flatten-doc.xsl \
	-xml \
	-out $@

$(FLATTENED_DIR)/compartment.rng : | $(FLATTENED_DIR)

$(FLATTENED_DIR)/compartment.rng : $(DOC_SCHEMA_SOURCE_DIR)/compartment.rng $(CPT_SCMTA)
	java org.apache.xalan.xslt.Process \
	-in $(DOC_SCHEMA_SOURCE_DIR)/compartment.rng \
	-xsl $(XSL)/flatten-doc.xsl \
	-xml \
	-out $@

$(FLATTENED_DIR)/rk4tau.rng : | $(FLATTENED_DIR)

$(FLATTENED_DIR)/rk4tau.rng : $(DOC_SCHEMA_SOURCE_DIR)/rk4tau.rng
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/flatten-doc.xsl \
	-xml \
	-out $@

$(FLATTENED_DIR)/rk4ode.rng : | $(FLATTENED_DIR)

$(FLATTENED_DIR)/rk4ode.rng : $(DOC_SCHEMA_SOURCE_DIR)/rk4ode.rng
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/flatten-doc.xsl \
	-xml \
	-out $@

$(FLATTENED_DIR)/odie.rng : | $(FLATTENED_DIR)

$(FLATTENED_DIR)/odie.rng : $(DOC_SCHEMA_SOURCE_DIR)/odie.rng
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/flatten-doc.xsl \
	-xml \
	-out $@

PREEN_LIST += $(DOC_SCHEMA_SOURCE_DIR)/*.bak
