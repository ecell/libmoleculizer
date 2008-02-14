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

$(VAR_SUITE_OUTPUT_DIR) : $(OUTPUT_DIR)
	mkdir $@

# This target exists so that all the directories can be created before the
# perl script that generates the substituted moleculizer input files runs.
var-suite-dirs : | $(VAR_SUITE_OUTPUT_DIR)

TEMPLATE_INPUT := $(VAR_SUITE_OUTPUT_DIR)/moleculizer-in.xml

$(VAR_SUITE_OUTPUT_DIR)/variations-setup : VARS_FILE := $(VARIATIONS_FILE)
$(VAR_SUITE_OUTPUT_DIR)/variations-setup : TPLT_INPUT := $(TEMPLATE_INPUT)
$(VAR_SUITE_OUTPUT_DIR)/variations-setup : RT_DIR := $(VAR_SUITE_OUTPUT_DIR)
$(VAR_SUITE_OUTPUT_DIR)/variations-setup : feedbacks-setup \
	$(TEMPLATE_INPUT) \
	$(VARIATIONS_FILE) \
	$(VAR_SUITE_OUTPUT_DIR)/dose-response-dirs
	$(SCRIPT)/do-substitutions.pl $(TPLT_INPUT) \
	$(VARS_FILE) $(RT_DIR) moleculizer-in.xml

include $(TMP)/variations.mk

TDDR_DIR := $(VAR_SUITE_OUTPUT_DIR)/tddrs

$(TDDR_DIR) : | $(VAR_SUITE_OUTPUT_DIR)
	mkdir $@

include $(TMP)/summary-responses.mk