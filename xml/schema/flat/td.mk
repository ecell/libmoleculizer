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

DOT := $(DOT)/flat

# These flattened schema-docs are the targets of this makefile.
# Here, "flattening" means satisfying includes, not refs.
INVENTORY := $(DOT)/moleculizer.rng \
	$(DOT)/odie.rng \
	$(DOT)/rk4tau.rng \
	$(DOT)/rk4ode.rng \
	$(DOT)/mzrState.rng

# Schema files included in moleculizer.rng and mzrState.rng.
MOLECULIZER_SCHEMATA := dimer.rng \
	gpa.rng \
	modKinase.rng \
	mol.rng \
	moleculizer.rng \
	nucEx.rng \
	plex.rng \
	scaffold.rng \
	stoch.rng \
	bndKinase.rng

# The locations of the MOLECULIZER_SCHEMATA.
MZR_SCMTA := $(addprefix $(SCHEMA)/schema-doc/,$(MOLECULIZER_SCHEMATA))

$(DOT)/target : $(INVENTORY)

$(DOT)/moleculizer.rng : $(SCHEMA)/schema-doc/moleculizer.rng $(MZR_SCMTA)
	java org.apache.xalan.xslt.Process \
	-in $(SCHEMA)/schema-doc/moleculizer.rng \
	-xsl $(XSL)/flatten-doc.xsl \
	-xml \
	-out $@

$(DOT)/mzrState.rng : $(SCHEMA)/schema-doc/mzrState.rng $(MZR_SCMTA)
	java org.apache.xalan.xslt.Process \
	-in $(SCHEMA)/schema-doc/mzrState.rng \
	-xsl $(XSL)/flatten-doc.xsl \
	-xml \
	-out $@

$(DOT)/rk4tau.rng : $(SCHEMA)/schema-doc/rk4tau.rng
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/flatten-doc.xsl \
	-xml \
	-out $@

$(DOT)/rk4ode.rng : $(SCHEMA)/schema-doc/rk4ode.rng
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/flatten-doc.xsl \
	-xml \
	-out $@

$(DOT)/odie.rng : $(SCHEMA)/schema-doc/odie.rng
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/flatten-doc.xsl \
	-xml \
	-out $@

CLEAN_LIST := $(CLEAN_LIST) $(INVENTORY)

PREEN_LIST := $(PREEN_LIST) $(DOT)/*~ $(DOT)/*.bak

DOT := $(call dotdot,$(DOT))
