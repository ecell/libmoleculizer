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

VERSION_NAME := version7

VERSION_DIR := $(DEMO_DIR)/$(VERSION_NAME)

INPUT_FILE := $(VERSION_DIR)/$(DEMO_NAME).xml

DOC_INPUT_HTML := $(VERSION_DIR)/source.html

OUTPUT_DIR := $(VERSION_DIR)/$(DEMO_NAME).out

$(VERSION_DIR)/target : $(OUTPUT_DIR)/plots-done $(DOC_INPUT_HTML)

$(OUTPUT_DIR) :
	mkdir $@

# Run the simulation in the output directory, generating many output files.
$(OUTPUT_DIR)/simulation-done : $(INPUT_FILE) | $(OUTPUT_DIR)
	echo "Started at" `date` > $@
	cp $< $(@D)
	cd $(@D) \
	&& cptzr < $(<F)
	echo "Finished at" `date` >> $@

# Makes time detail of five-total-membranes.dmp.
$(OUTPUT_DIR)/five-total-membranes-sec15-20.dmp : $(OUTPUT_DIR)/simulation-done
	head -n 1 $(@D)/five-total-membranes.dmp > $@
	tail -n 1700 $(@D)/five-total-membranes.dmp | head -n 500 >> $@

# Makes compartment details of above.
$(OUTPUT_DIR)/five-total-membranes-sec15-20-c2.dmp : $(OUTPUT_DIR)/five-total-membranes-sec15-20.dmp
	cut -f1,2 $< > $@

$(OUTPUT_DIR)/five-total-membranes-sec15-20-c7.dmp : $(OUTPUT_DIR)/five-total-membranes-sec15-20.dmp
	cut -f1,7 $? > $@

# Fixes column headers.  This needs to be dealt with in Moleculizer's code,
# one needs a way of optionally giving special column headers in the
# dump-stream specification.
#
# Also, putting in a back up step here, since the line numbers can change
# and it can be necessary to re-run this target to get it right.
$(OUTPUT_DIR)/headers-fixed : \
	$(OUTPUT_DIR)/five-total.dmp.bak \
	$(OUTPUT_DIR)/four-five-total.dmp.bak
	echo "#sim-time	five-total	five-total:membrane	five-total:cytoplasm" > $(@D)/tmp
	tail -n 3200 $(@D)/five-total.dmp >> $(@D)/tmp
	mv $(@D)/tmp $(@D)/five-total.dmp
	echo "#sim-time	four-five	four-five:membrane	four-five:cytoplasm" > $(@D)/tmp
	tail -n 3200 $(@D)/four-five-total.dmp >> $(@D)/tmp
	mv $(@D)/tmp $(@D)/four-five-total.dmp
	touch $@

# Backup because line numbers above can go bad.
$(OUTPUT_DIR)/five-total.dmp.bak : $(OUTPUT_DIR)/simulation-done
	cp $(@D)/five-total.dmp $@

# Backup because line numbers above can go bad.
$(OUTPUT_DIR)/four-five-total.dmp.bak : $(OUTPUT_DIR)/simulation-done
	cp $(@D)/four-five-total.dmp $@

# Plot all the .dmp files.
$(OUTPUT_DIR)/plots-done : \
	$(OUTPUT_DIR)/simulation-done \
	$(OUTPUT_DIR)/five-total-membranes-sec15-20.dmp \
	$(OUTPUT_DIR)/five-total-membranes-sec15-20-c2.dmp \
	$(OUTPUT_DIR)/five-total-membranes-sec15-20-c7.dmp \
	$(OUTPUT_DIR)/headers-fixed
	cd $(@D) && plot-dmp-files
	touch $@

# Make documented html version of input file.
# 
# I'm somewhat uncomfortable using this relative path for the documentation
# URL.
$(DOC_INPUT_HTML) : FLAT_DOC_URL := ../../../doc/static-doc/
$(DOC_INPUT_HTML) : $(INPUT_FILE)
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/arg-any2static-doc.xsl \
	-param flat-doc-url $(FLAT_DOC_URL) \
	-param caption $(<F) \
	-html \
	-out $@

$(VERSION_DIR)/clean : CLN := $(OUTPUT_DIR) $(DOC_INPUT_HTML)
$(VERSION_DIR)/clean :
	rm -rf $(CLN)

CLEAN_LIST += $(OUTPUT_DIR) $(DOC_INPUT_HTML)

PREEN_LIST += $(VERSION_DIR)/*~ $(VERSION_DIR)/*.bak

