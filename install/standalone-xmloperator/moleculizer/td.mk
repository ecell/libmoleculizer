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

DOT := $(DOT)/moleculizer

# Build-time target: schemas and xsl transformers.
$(DOT)/target : $(DOT)/schema/target $(DOT)/xsl/target 

# Install-time target: user configuration files for user.xmloperator.
$(DOT)/install-target : $(DOT)/mzr-defaults.xml \
	$(DOT)/xmlo_rc \
	$(DOT)/user.molrc

$(DOT)/user.molrc : $(INSTALL)/misc/user.molrc
	cp $< $@

$(DOT)/xmlo_rc : $(INSTALL)/misc/xmlo_rc
	cp $< $@

$(DOT)/mzr-defaults.xml : $(INSTALL)/misc/mzr-defaults.xml
	cp $< $@

INSTALL_CLEAN_LIST := $(INSTALL_CLEAN_LIST) \
	$(DOT)/user.molrc \
	$(DOT)/xmlo_rc \
	$(DOT)/mzr-defaults.xml \

PREEN_LIST := $(PREEN_LIST) $(DOT)/*~ $(DOT)/*.bak

include $(DOT)/schema/td.mk $(DOT)/xsl/td.mk

DOT := $(call dotdot,$(DOT))
