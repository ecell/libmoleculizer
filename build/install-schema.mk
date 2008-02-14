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

# Rules to put schema files, ultimately destined to go into user.xmloperator,
# into the staging area for user installation material.  This staging area
# is included in the binary release.

# The directory where schema files go.
INS_SCHEMA_DIR := $(INS_TGT_DIR)/schema

# This directory is copied from install/xmloperator when user.xmloperator is
# copied, so it doesn't have to be mkdir'ed.
$(INS_SCHEMA_DIR) : | $(INS_TGT_DIR)

# The schemas that are generated from doc-schemas.  DOC_SCHEMAS comes from
# schema .mk.
INS_DOC_SCHEMAS := $(addprefix $(INS_SCHEMA_DIR)/,$(DOC_SCHEMAS))

$(INS_DOC_SCHEMAS) : | $(INS_SCHEMA_DIR)

$(INS_DOC_SCHEMAS) : $(INS_SCHEMA_DIR)/% : $(DOC_REMOVED_DIR)/%
	cp $< $@

# The schemas that are just written by hand.  FIXED_SCHEMAS comes from
# schema.mk.
INS_FIXED_SCHEMAS := $(addprefix $(INS_SCHEMA_DIR)/,$(FIXED_SCHEMAS))

$(INS_FIXED_SCHEMAS) : | $(INS_SCHEMA_DIR)

# Right now, FIXED_SCHEMA_SOURCE_DIR is in the source area; it needs to be
# moved to the top-level target directory for xml, and the need to be copied
# there from the source area, in order for this to work in the binary release.
$(INS_FIXED_SCHEMAS) : $(INS_SCHEMA_DIR)/% : $(FIXED_SCHEMA_SOURCE_DIR)/%
	cp $< $@

# Add the schema targets to the overall target.
$(INS_XMLOPERATOR)/.configurable : $(INS_FIXED_SCHEMAS) \
	$(INS_DOC_SCHEMAS)
