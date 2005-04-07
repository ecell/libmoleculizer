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

DOT := $(DOT)/schema

# These are the schemas that are written in rngdoc form and have to
# be converted for use by xmlOperator.
#
# The rngdoc.rng schema for the documented form of the rng schema is not
# itself transformed.  (xmlOperator seemed to have some sort of difficulty
# with the output of that operation.)
#
# May want to remove rk4ode.rng, since that app may not go into the release.
# Though users won't want to edit moleculizer-state, it will still need to
# be a document type in xmloperator.
DOC_SCHEMA := $(DOT)/dimer.rng \
	$(DOT)/gpa.rng \
	$(DOT)/modKinase.rng \
	$(DOT)/bndKinase.rng \
	$(DOT)/ftr.rng \
	$(DOT)/mol.rng \
	$(DOT)/moleculizer.rng \
	$(DOT)/mzrState.rng \
	$(DOT)/nucEx.rng \
	$(DOT)/odie.rng \
	$(DOT)/plex.rng \
	$(DOT)/rk4ode.rng \
	$(DOT)/rk4tau.rng \
	$(DOT)/scaffold.rng \
	$(DOT)/stoch.rng

# These are schema that are not generated from rngdoc.  Leaving out
# rngdoc.rng, since causal users won't need to edit those.  This is needed
# so that a user could gracefully change her server, for example, but it could
# really be omitted, too.
NO_DOC_SCHEMA := $(DOT)/mzr-defaults.rng \
	$(DOT)/rngdoc.rng

$(DOT)/target : $(DOC_SCHEMA) $(NO_DOC_SCHEMA)

# Make the local schema files by removing documentation from the original
# documented schema files.
$(DOC_SCHEMA) : $(DOT)/% : $(SCHEMA)/schema-doc/%
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/doc2rng.xsl \
	-xml \
	-out $@

# Make local copies of schema files that are note produced from documented
# schema files.
$(NO_DOC_SCHEMA) : $(DOT)/% : $(SCHEMA)/fixed/%
	cp $< $@

CLEAN_LIST := $(CLEAN_LIST) $(DOC_SCHEMA) $(NO_DOC_SCHEMA)

PREEN_LIST := $(PREEN_LIST) $(DOT)/*~ $(DOT)/*.bak

DOT := $(call dotdot,$(DOT))
