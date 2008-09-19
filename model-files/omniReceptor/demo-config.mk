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

DEMO_NAME := omniReceptor

DEMO_DIR := $(DEMO)/$(DEMO_NAME)

DEMO_INPUT := $(DEMO_DIR)/$(DEMO_NAME).xml

# Output directory for the usual simulation run.
OUTPUT_DIR := $(DEMO_DIR)/$(DEMO_NAME).out

# Output directories for runs with look-ahead turned on.
DEPTH_ONE_DIR := $(DEMO_DIR)/depth-one.out
DEPTH_TWO_DIR := $(DEMO_DIR)/depth-two.out

# The default target, which just runs the simulation in the usual way.
$(DEMO_DIR)/target : $(OUTPUT_DIR)/simulation-done

# Tests the lookahead feature, running the simulation with one and two levels
# of lookahead.
$(DEMO_DIR)/lookahead : $(DEPTH_ONE_DIR)/simulation-done \
	$(DEPTH_TWO_DIR)/simulation-done

$(OUTPUT_DIR) :
	mkdir $@

# Simulation done without species look-ahead; an ordinary moleculizer
# simulation.
$(OUTPUT_DIR)/simulation-done : $(DEMO_INPUT) | $(OUTPUT_DIR)
	echo "Started at" `date` > $@
	cp $< $(@D)
	cd $(@D) \
	&& moleculizer < $(<F) \
	&& plot-dmp-files
	echo "Finished at" `date` >>$@

$(DEPTH_ONE_DIR) :
	mkdir $@

# Simulation done with 1 generation look-ahead.
$(DEPTH_ONE_DIR)/simulation-done : $(DEMO_INPUT) | $(DEPTH_ONE_DIR) 
	echo "Started at" `date` > $@
	cp $< $(@D)
	cd $(@D) \
	&& moleculizer -d 1 < $(<F) \
	&& plot-dmp-files
	echo "Finished at" `date` >>$@

$(DEPTH_TWO_DIR) :
	mkdir $@

# Simulation done with 2 generations look-ahead.
$(DEPTH_TWO_DIR)/simulation-done : $(DEMO_INPUT) | $(DEPTH_TWO_DIR) 
	echo "Started at" `date` > $@
	cp $< $(@D)
	cd $(@D) \
	&& moleculizer -d 2 < $(<F) \
	&& plot-dmp-files
	echo "Finished at" `date` >>$@

# For redoing the demo repeatedly.
$(DEMO_DIR)/clean : CLN := $(OUTPUT_DIR) $(DEPTH_ONE_DIR) $(DEPTH_TWO_DIR)
$(DEMO_DIR)/clean :
	rm -rf $(CLN)

CLEAN_LIST +=  $(OUTPUT_DIR) $(DEPTH_ONE_DIR) $(DEPTH_TWO_DIR)

PREEN_LIST += $(DEMO_DIR)/*~ $(DEMO_DIR)/*.bak
