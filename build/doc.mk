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

# Target directory for documentation.
$(DOC) :
	mkdir $@

##############

# Target directory for doxygen-generated source code documentation.
DOXYGEN_TARGET_DIR := $(DOC)/doxygen

# Note that this configuration file must be updated when things move around.
DOXYGEN_CONFIGURATION_FILE := $(DOCS)/doxygen-config

$(DOXYGEN_TARGET_DIR) : | $(DOC)
	mkdir $@

$(DOXYGEN_TARGET_DIR)/.target : | $(DOXYGEN_TARGET_DIR)
	doxygen $(DOXYGEN_CONFIGURATION_FILE)
	touch $@

# The general target for doxygen source code documentation.
src-doc : $(DOXYGEN_TARGET_DIR)/.target

##############

# A few miscellaneous html and pdf doc files.

$(DOC)/index.html : $(DOCS)/index.html | $(DOC)
	cp $< $@

$(DOC)/mzr_bio.pdf : $(DOCS)/mzr_bio.pdf | $(DOC)
	cp $< $@

usr-doc : $(DOC)/index.html \
	$(DOC)/mzr_bio.pdf \

##############

# Include rules to generate HTML trees of documentation from schema-doc
# schemas.
include $(BUILD)/doc-static.mk

# Include rules to copy hand-generated html documentation to doc.
include $(BUILD)/doc-overflow.mk

##############

PREEN_LIST += $(DOCS)/*~ \
	$(DOCS)/*.bak \
