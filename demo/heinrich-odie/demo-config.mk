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

DEMO_NAME := heinrich-odie

DEMO_DIR := $(DEMO)/$(DEMO_NAME)

ODIE_INPUT_FILE := $(DEMO_DIR)/$(DEMO_NAME).xml

# In this case, html translation of the odie input file.
DOC_INPUT_HTML := $(DEMO_DIR)/source.html

ODIE_OUTPUT_DIR := $(DEMO_DIR)/$(DEMO_NAME).out

MZR_INPUT_FILE := $(DEMO_DIR)/heinrich.xml

MZR_OUTPUT_DIR := $(DEMO_DIR)/heinrich.out

# We convert the state dump to odie input, for comparison with the odie input
# that is actually used.  But we retain the edited odie input, so that the
# hand editing doesn't have to be done every time.
$(DEMO_DIR)/target : $(ODIE_OUTPUT_DIR)/simulation-done \
	$(DEMO_DIR)/heinrich-odie-rough.xml \
	$(DOC_INPUT_HTML)

$(ODIE_OUTPUT_DIR) :
	mkdir $@

$(MZR_OUTPUT_DIR) :
	mkdir $@

# Run the odie simulation, which comes via XSLT and hand editing from the
# moleculizer state dump.
$(ODIE_OUTPUT_DIR)/simulation-done : $(ODIE_INPUT_FILE) | $(ODIE_OUTPUT_DIR)
	echo "Started at " `date` > $@
	cp $< $(@D)
	cd $(@D) \
	&& odie < $(<F) \
	&& plot-dmp-files-as-cnc
	echo "Finished at " `date` >> $@

# Run the moleculizer simulation, which dumps state.
# 
# I don't know what this parametrizer run is all about...
$(MZR_OUTPUT_DIR)/simulation-done : $(MZR_INPUT_FILE) | $(MZR_OUTPUT_DIR)
	echo "Started at " `date` > $@
	cp $< $(@D)
	cd $(@D) \
	&& moleculizer < $(<F) \
	&& parametrizer heinrich-state.xml < $(<F) > param-dump.xml \
	&& plot-dmp-files
	echo "Finished at " `date` >> $@

# Convert the state dump into odie input.  The resulting odie input needs a
# little hand-editing to set some odie-specific parameters.
$(DEMO_DIR)/heinrich-odie-rough.xml : $(MZR_OUTPUT_DIR)/simulation-done
	state2odie $(<D)/param-dump.xml $@

# Make documented html version of input file.
# 
# I'm somewhat uncomfortable using this relative path for the documentation
# URL.
$(DOC_INPUT_HTML) : FLAT_DOC_URL := ../../doc/static-doc/
$(DOC_INPUT_HTML) : $(ODIE_INPUT_FILE)
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/arg-any2static-doc.xsl \
	-param flat-doc-url $(FLAT_DOC_URL) \
	-param caption $(<F) \
	-html \
	-out $@

# For redoing individual demo several times.
$(DEMO_DIR)/clean : CLN := $(ODIE_OUTPUT_DIR) \
	$(MZR_OUTPUT_DIR) \
	$(DEMO_DIR)/heinrich-odie-rough.xml \
	$(DOC_INPUT_HTML)
$(DEMO_DIR)/clean :
	rm -rf $(CLN)

CLEAN_LIST += $(ODIE_OUTPUT_DIR) \
	$(MZR_OUTPUT_DIR) \
	$(DEMO_DIR)/heinrich-odie-rough.xml \
	$(DOC_INPUT_HTML)

PREEN_LIST += $(DEMO_DIR)/*~ $(DEMO_DIR)/*.bak
