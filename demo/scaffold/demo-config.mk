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

DEMO_NAME := scaffold

DEMO_DIR := $(DEMO)/$(DEMO_NAME)

# This input file is used in the "with alpha" simulation.
# Another input file, which must be created by hand from this one,
# is used in the basal response simulation.
INPUT_FILE := $(DEMO_DIR)/$(DEMO_NAME).xml

DOC_INPUT_HTML := $(DEMO_DIR)/source.html

OUTPUT_DIR := $(DEMO_DIR)/$(DEMO_NAME).out

# The simulation with added alpha factor.
WITH_ALPHA_DIR := $(OUTPUT_DIR)/with.out

# The basal simulation, without added alpha factor.
WITHOUT_ALPHA_DIR := $(OUTPUT_DIR)/without.out

# The input file for the basal simulation.  This must be generated
# by hand from the input file for the ordinary simulation, so it
# is not deleted by clean.
WITHOUT_INPUT_FILE := $(DEMO_DIR)/without-alpha.xml

# Standard target; doesn't do the basal simulation.
$(DEMO_DIR)/target : $(DEMO_DIR)/with $(DOC_INPUT_HTML)

# Do both standard target and basal simulation.
$(DEMO_DIR)/both : $(DEMO_DIR)/target $(DEMO_DIR)/without

# Do simulation with alpha factor added.
$(DEMO_DIR)/with : $(WITH_ALPHA_DIR)/simulation-done

# Do basal simulation without alpha factor added.
$(DEMO_DIR)/without : $(WITHOUT_ALPHA_DIR)/simulation-done

# Make the overall output directory.
$(OUTPUT_DIR) :
	mkdir $@

# Make the output directory for the normal simulation with alpha factor.
$(WITH_ALPHA_DIR) : $(OUTPUT_DIR)
	mkdir $@

# Copy the input file for the normal simulation to its output directory and
# run the simulation.
$(WITH_ALPHA_DIR)/simulation-done : $(INPUT_FILE) | $(WITH_ALPHA_DIR)
	echo "Started at" `date` > $@
	cp $< $(@D)
	cd $(@D) \
	&& moleculizer < $(<F) \
	&& plot-dmp-files
	echo "Finished at" `date` >> $@

# Rule to remind the user that she needs to delete the addition of
# alpha factor from the simulation input by hand.
$(WITHOUT_INPUT_FILE) : $(INPUT_FILE)

# Make the output directory for the basal simulation without alpha factor.
$(WITHOUT_ALPHA_DIR) : $(OUTPUT_DIR)
	mkdir $@

# Copy the input file for the basal simulation to its output directory and
# run the simulation.
$(WITHOUT_ALPHA_DIR)/simulation-done : $(WITHOUT_INPUT_FILE) | $(WITHOUT_ALPHA_DIR)
	echo "Started at" `date` > $@
	cp $< $(@D)
	cd $(@D) \
	&& moleculizer < $(<F) \
	&& plot-dmp-files
	echo "Finished at" `date` >> $@

# Translate the input file to html.
$(DOC_INPUT_HTML) : FLAT_DOC_URL := ../../doc/static-doc/
$(DOC_INPUT_HTML) : $(INPUT_FILE)
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/arg-mzr2static-doc.xsl \
	-param flat-doc-url $(FLAT_DOC_URL) \
	-param caption $(<F) \
	-html \
	-out $@

# To assist running the demo several times.
$(DEMO_DIR)/clean : CLN := $(OUTPUT_DIR) $(DOC_INPUT_HTML)
$(DEMO_DIR)/clean :
	rm -rf $(CLN)

CLEAN_LIST += $(OUTPUT_DIR) $(DOC_INPUT_HTML)

PREEN_LIST += $(DEMO_DIR)/*~ $(DEMO_DIR)/*.bak
