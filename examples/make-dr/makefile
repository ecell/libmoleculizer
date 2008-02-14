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

include build/navigation.mk

target : $(OUTPUT_DIR)/feedbacks

include build/build.mk

include build/script.mk

$(OUTPUT_DIR) :
	mkdir $@

TEMPLATE_INPUT := ./moleculizer-in.xml

$(VARIATIONS_FILE) : VARIATION_COUNT := 5
$(VARIATIONS_FILE) : VARIATION_FRACTION := 0.3
$(VARIATIONS_FILE) : moleculizer-pops
	nu -load $(SCRIPT)/randomize-pops.scm -input $< \
	-output $@ -count $(VARIATION_COUNT) -fraction $(VARIATION_FRACTION)

feedbacks-setup : TPLT_INPUT := $(TEMPLATE_INPUT)
feedbacks-setup : FBKS_FILE := $(FEEDBACKS_FILE)
feedbacks-setup : $(TEMPLATE_INPUT) \
	$(FEEDBACKS_FILE) \
	var-suite-dirs 
	$(SCRIPT)/do-substitutions.pl $(TPLT_INPUT) \
	$(FBKS_FILE) $(OUTPUT_DIR) moleculizer-in.xml

include $(TMP)/feedbacks.mk

clean :
	rm -rf $(OUTPUT_DIR)
	rm -f $(VARIATIONS_FILE)

preen : clean
	rm -f *~
	rm -f *.bak
