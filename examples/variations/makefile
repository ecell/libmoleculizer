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

target : variations-target

test-target :

include build/navigation.mk

include build/build.mk

include build/script.mk

$(VARIATIONS_FILE) : VARIATION_COUNT := 10
$(VARIATIONS_FILE) : VARIATION_FRACTION := 0.3
$(VARIATIONS_FILE) : moleculizer-pops
	nu -load $(SCRIPT)/randomize-pops.scm -input $< \
	-output $@ -count $(VARIATION_COUNT) -fraction $(VARIATION_FRACTION)

variations-target : $(VARIATIONS_FILE)

TEMPLATE_INPUT := ./moleculizer-in.xml

variations-setup : VARS_FILE := $(VARIATIONS_FILE)
variations-setup : TPLT_INPUT := $(TEMPLATE_INPUT)
variations-setup : $(TEMPLATE_INPUT) $(VARIATIONS_FILE) \
	dose-response-dirs
	$(SCRIPT)/do-substitutions.pl $(TPLT_INPUT) \
	$(VARS_FILE) $(OUTPUT_DIR) moleculizer-in.xml

$(OUTPUT_DIR) :
	mkdir $@

clean : VARS := $(VARIATIONS_FILE)
clean :
	rm -rf $(OUTPUT_DIR)
	rm -f $(VARS)

preen : clean
	rm -f *~
	rm -f *.bak

include $(TMP)/variations.mk

TDDR_DIR := $(OUTPUT_DIR)/tddrs

$(TDDR_DIR) : | $(OUTPUT_DIR)
	mkdir $@

include $(TMP)/summary-responses.mk