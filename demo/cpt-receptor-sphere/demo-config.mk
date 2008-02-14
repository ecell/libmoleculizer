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

DEMO_NAME := cpt-receptor-sphere

DEMO_DIR := $(DEMO)/$(DEMO_NAME)

VERSIONS := version1 \
	version2 \
	version3 \
	version4 \
	version5 \
	version6 \
	version7

VERSION_TARGETS := $(addprefix $(DEMO_DIR)/,$(addsuffix /target,$(VERSIONS)))
VERSION_CLEANS := $(addprefix $(DEMO_DIR)/,$(addsuffix /clean,$(VERSIONS)))

# Run the last version of the demo.
$(DEMO_DIR)/target : $(DEMO_DIR)/version7/target

# Run all versions of the demo.
$(DEMO_DIR)/all : $(VERSION_TARGETS)

# Clean all versions.
$(DEMO_DIR)/clean : $(VERSION_CLEANS)

PREEN_LIST += $(DEMO_DIR)/*~ $(DEMO_DIR)/*.bak

# Configuration files for all the versions. They differ slightly from one
# another.
VERSION_CONFIGS := $(addprefix $(DEMO_DIR)/,$(addsuffix /demo-config.mk,$(VERSIONS)))

include $(VERSION_CONFIGS)
