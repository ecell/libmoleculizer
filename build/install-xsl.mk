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

# Rules to copy xslt transformations, ultimately destined to go into
# user.xmloperator, into the staging area for user installation material.
# This staging area is included in the binary release.

# XSLT transformations used by xmloperator to generate live documentation.
XMLO_XSLT := moleculizer2doc.xsl \
	mzr-tests.xsl \
	mol-tests.xsl \
	plex-tests.xsl \
	dimer-tests.xsl \
	gpa-tests.xsl \
	nuc-ex-tests.xsl \
	kinase-tests.xsl \
	scaffold-tests.xsl \
	stateDump2doc.xsl \
	moleculizer2static-doc.xsl

INS_XSL_DIR := $(INS_TGT_DIR)/xsl

INS_XSL_FILES := $(addprefix $(INS_XSL_DIR)/,$(XMLO_XSLT))

$(INS_XSL_FILES) : $(INS_XSL_DIR)/% : $(XSL_TARGET_DIR)/%
	cp $< $@

# Add the xslt targets to the overall target.
$(INS_XMLOPERATOR)/.configurable : $(INS_XSL_FILES)


