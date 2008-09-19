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

DEMO_NAME := simple-rk4tau

DEMO_DIR := $(DEMO)/$(DEMO_NAME)

MZR_INPUT_FILE := $(DEMO_DIR)/simple.xml

MZR_OUTPUT_DIR := $(DEMO_DIR)/simple.out

RK4TAU_TMP_INPUT_FILE := $(DEMO_DIR)/$(DEMO_NAME)-tmp.xml

RK4TAU_INPUT_FILE := $(DEMO_DIR)/$(DEMO_NAME).xml

RK4TAU_OUTPUT_DIR := $(DEMO_DIR)/$(DEMO_NAME).out

$(DEMO_DIR)/target : $(RK4TAU_OUTPUT_DIR)/simulation-done

$(DEMO_DIR)/input : $(RK4TAU_TMP_INPUT_FILE)

$(RK4TAU_OUTPUT_DIR) :
	mkdir $@

# Run the rk4tau simulation, using reaction network generated from moleculizer
# state dump.
$(RK4TAU_OUTPUT_DIR)/simulation-done : $(RK4TAU_INPUT_FILE) | $(RK4TAU_OUTPUT_DIR)
	echo "Started at" `date` > $@
	cp $< $(@D)
	cd $(@D) \
	&& rk4tau < $(<F) \
	&& plot-dmp-files
	echo "Finished at" `date` >> $@

# The rk4tau input file generated from the moleculizer state dump requires
# a little hand editing on rk4tau-specific parameters.  With this approach,
# one doesn't have to do the hand-editing every time.
$(RK4TAU_TMP_INPUT_FILE) : $(MZR_OUTPUT_DIR)/simulation-done
	state2rk4tau $(<D)/param-dump.xml $@

$(MZR_OUTPUT_DIR) :
	mkdir $@

# Run the moleculizer simulation, generating state dump that is translated
# into rk4tau input.
$(MZR_OUTPUT_DIR)/simulation-done : $(MZR_INPUT_FILE) | $(MZR_OUTPUT_DIR)
	echo "Started at" `date` > $@
	cp $< $(@D)
	cd $(@D) \
	&& moleculizer < $(<F) \
	&& parametrizer state-dump.xml < $(<F) > param-dump.xml \
	&& plot-dmp-files
	echo "Finished at" `date` >> $@

# For running the demo repeatedly.
$(DEMO_DIR)/clean : CLN := $(RK4TAU_OUTPUT_DIR) \
	$(MZR_OUTPUT_DIR) \
	$(RK4TAU_TMP_INPUT_FILE)
$(DEMO_DIR)/clean :
	rm -rf $(CLN)

CLEAN_LIST += $(RK4TAU_OUTPUT_DIR) \
	$(MZR_OUTPUT_DIR) \
	$(RK4TAU_TMP_INPUT_FILE)

PREEN_LIST += $(DEMO_DIR)/*~ $(DEMO_DIR)/*.bak
