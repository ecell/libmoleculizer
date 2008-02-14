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

XSL_FILES := any2doc.xsl \
	any2static-doc.xsl \
	arg-mzr2static-doc.xsl \
	dimer-tests.xsl \
	doc2rng.xsl \
	dump2odie.xsl \
	dump2sbml.xsl \
	dump2tau.xsl \
	editAttr.xsl \
	edit-value.xsl \
	edit-xpth.xsl \
	flatten-doc.xsl \
	gpa-tests.xsl \
	inv2html.xsl \
	kinase-tests.xsl \
	moleculizer2doc.xsl \
	moleculizer2static-doc.xsl \
	mol-tests.xsl \
	mzr2form.xsl \
	mzr2range.xsl \
	mzr-tests.xsl \
	nuc-ex-tests.xsl \
	odie2static-doc.xsl \
	origPops.xsl \
	plex-tests.xsl \
	rk4tau2static-doc.xsl \
	rng2doc.xsl \
	rng-doc-ndx.xsl \
	rng-doc-toc.xsl \
	rngpth.xsl \
	scaffold-tests.xsl \
	state2static-doc.xsl \
	stateDump2doc.xsl \
	static-doc-ndx.xsl \
	static-doc-toc.xsl \
	static-doc.xsl

XSL_TARGET_DIR := $(XML)/xsl

$(XSL_TARGET_DIR) : | $(XML)
	mkdir $@

XSL_TARGET_FILES := $(addprefix $(XSL_TARGET_DIR)/,$(XSL_FILES))

$(XSL_TARGET_FILES) : $(XSL_TARGET_DIR)/% : $(XSL)/% | $(XSL_TARGET_DIR)
	cp $< $@

xml-target : $(XSL_TARGET_FILES)

PREEN_LIST += $(XSL)/*.bak