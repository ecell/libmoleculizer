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

DEMO_NAME := heinrich-rk4tau

DEMO_DIR := $(DEMO)/$(DEMO_NAME)

RK4TAU_INPUT_FILE := $(DEMO_DIR)/$(DEMO_NAME).xml

# In this case, the html translation of the rk4tau input.
DOC_INPUT_HTML := $(DEMO_DIR)/source.html

RK4TAU_OUTPUT_DIR := $(DEMO_DIR)/$(DEMO_NAME).out

MZR_INPUT_FILE := $(DEMO_DIR)/heinrich.xml

MZR_OUTPUT_DIR := $(DEMO_DIR)/heinrich.out

$(DEMO_DIR)/target : $(RK4TAU_OUTPUT_DIR)/simulation-done \
	$(DEMO_DIR)/heinrich-rk4tau-rough.xml \
	$(DOC_INPUT_HTML)

$(RK4TAU_OUTPUT_DIR) :
	mkdir $@

# Run the rk4tau simulation, which comes via XSLT and hand editing from the
# moleculizer state dump.
$(RK4TAU_OUTPUT_DIR)/simulation-done : $(RK4TAU_INPUT_FILE) | $(RK4TAU_OUTPUT_DIR)
	echo "Started at" `date` > $@
	cp $< $(@D)
	cd $(@D) \
	&& rk4tau < $(<F) \
	&& plot-dmp-files-as-cnc
	echo "Finished at" `date` >> $@

$(MZR_OUTPUT_DIR) :
	mkdir $@

# Run the moleculizer simulation, wich dumps state.  I don't see why this
# parametrizer run is here.
$(MZR_OUTPUT_DIR)/simulation-done : $(MZR_INPUT_FILE) | $(MZR_OUTPUT_DIR)
	echo "Started at" `date` > $@
	cp $< $(@D)
	cd $(@D) \
	&& moleculizer < $(<F) \
	&& parametrizer heinrich-state.xml < $(<F) > param-dump.xml \
	&& plot-dmp-files
	echo "Finished at" `date` >> $@

# Convert moleculizer state dump into rk4tau input.  The resulting rk4tau
# input needs a little hand editing to set some rk4tau-specific parameters.
# 
# We do this conversion for comparison with the edited rk4tau input that is
# actually used.  But we retain the edited rk4tau input, so that hand editing
# doesn't ahve to be done every time.
$(DEMO_DIR)/heinrich-rk4tau-rough.xml : $(MZR_OUTPUT_DIR)/simulation-done
	state2rk4tau $(<D)/param-dump.xml $@

# Make documented html version of input file.
# 
# I'm somewhat uncomfortable using this relative path for the documentation
# URL.
$(DOC_INPUT_HTML) : FLAT_DOC_URL := ../../doc/static-doc/
$(DOC_INPUT_HTML) : $(RK4TAU_INPUT_FILE)
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/arg-any2static-doc.xsl \
	-param flat-doc-url $(FLAT_DOC_URL) \
	-param caption $(<F) \
	-html \
	-out $@

# For redoing individual demos several times.
$(DEMO_DIR)/clean : CLN := $(RK4TAU_OUTPUT_DIR) \
	$(MZR_OUTPUT_DIR) \
	$(DEMO_DIR)/heinrich-rk4tau-rough.xml \
	$(DOC_INPUT_HTML)
$(DEMO_DIR)/clean :
	rm -rf $(CLN)

CLEAN_LIST += $(RK4TAU_OUTPUT_DIR) \
	$(MZR_OUTPUT_DIR) \
	$(DEMO_DIR)/heinrich-rk4tau-rough.xml \
	$(DOC_INPUT_HTML)

PREEN_LIST += $(DEMO_DIR)/*~ $(DEMO_DIR)/*.bak
