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

# Get tools for hierarchical approach to make and
# principal landmarks in build tree.
include navigation.mk

# Get server configuration.
include server-config.mk

# Installation configurations: things that depend on where
# the moleculizer ball'o'wax is located (MOLECULIZER_DIR)
MOLECULIZER_DIR := $(shell pwd)
XMLO_DIR := $(MOLECULIZER_DIR)/xmloperator_2_3_0
FLAT_DOC_DIR := $(MOLECULIZER_DIR)/doc/static-doc
MISC_DIR := $(MOLECULIZER_DIR)/install/misc

# Generate the user.xmloperator file for users on the machine
# where moleculizer is installed.  This is also for developers
# who want to use their development area as an installation.
standalone-install : $(INSTALL)/standalone-xmloperator/install-target
	rm -rf user.xmloperator
	cp -a $(INSTALL)/standalone-xmloperator user.xmloperator

INSTALL_CLEAN_LIST := user.xmloperator
install-clean :
	rm -rf $(INSTALL_CLEAN_LIST)

other-things :

include $(INSTALL)/td.mk \
	demo/td.mk
