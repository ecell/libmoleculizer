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

DOT := $(DOT)/static-doc

# This absolute file path is needed to handle a special property of Xalan
# XSLT processing: paths used in the redirect feature are normally taken
# to be relative the file of the input file.  In these cases, that would
# be somewhere off in XML land.
#
# This absolute path doesn't actually go into any of the pages themselves,
# which are linked to each other by relative links.
# 
# CURDIR is defined automatically by make and works well with $(MAKE) -C.
TARGET_DIRECTORY = $(CURDIR)/$(DOT)

INVENTORY := $(DOT)/moleculizer-pages \
	$(DOT)/moleculizer-toc.html \
	$(DOT)/moleculizer-index.html \
	$(DOT)/moleculizer-state-pages \
	$(DOT)/moleculizer-state-toc.html \
	$(DOT)/moleculizer-state-index.html \
	$(DOT)/odie-pages \
	$(DOT)/odie-toc.html \
	$(DOT)/odie-index.html \
	$(DOT)/rk4tau-pages \
	$(DOT)/rk4tau-toc.html \
	$(DOT)/rk4tau-index.html

$(DOT)/target : $(INVENTORY)

# Not including dependency on the tranformer scripts here.

$(DOT)/moleculizer-pages : TD := $(TARGET_DIRECTORY)
$(DOT)/moleculizer-pages : $(SCHEMA)/flat/moleculizer.rng
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc.xsl \
	-param target-directory $(TD) \
	-param toc-doc-url moleculizer-toc.html \
	-param flat-doc-url $(DOC)/overflow-html \
	-param ndx-doc-url moleculizer-index.html \
	-html
	touch $@

$(DOT)/moleculizer-toc.html : $(SCHEMA)/flat/moleculizer.rng
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc-toc.xsl \
	-param caption "moleculizer-input table of contents" \
	-html > $@

$(DOT)/moleculizer-index.html : $(SCHEMA)/flat/moleculizer.rng
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc-ndx.xsl \
	-param caption "moleculizer-input element index" \
	-html > $@

$(DOT)/moleculizer-state-pages : TD := $(TARGET_DIRECTORY)
$(DOT)/moleculizer-state-pages : $(SCHEMA)/flat/mzrState.rng
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc.xsl \
	-param target-directory $(TD) \
	-param toc-doc-url moleculizer-state-toc.html \
	-param flat-doc-url $(DOC)/overflow-html \
	-param ndx-doc-url moleculizer-state-index.html \
	-html
	touch $@

$(DOT)/moleculizer-state-toc.html : $(SCHEMA)/flat/mzrState.rng
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc-toc.xsl \
	-param caption "moleculizer-state table of contents" \
	-html > $@

$(DOT)/moleculizer-state-index.html : $(SCHEMA)/flat/mzrState.rng
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc-ndx.xsl \
	-param caption "moleculizer-input element index" \
	-html > $@

$(DOT)/odie-pages : TD := $(TARGET_DIRECTORY)
$(DOT)/odie-pages : $(SCHEMA)/flat/odie.rng
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc.xsl \
	-param target-directory $(TD) \
	-param toc-doc-url odie-toc.html \
	-param flat-doc-url $(DOC)/overflow-html \
	-param ndx-doc-url odie-index.html \
	-html
	touch $@

$(DOT)/odie-toc.html : $(SCHEMA)/flat/odie.rng
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc-toc.xsl \
	-param caption "odie table of contents" \
	-html > $@

$(DOT)/odie-index.html : $(SCHEMA)/flat/odie.rng
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc-ndx.xsl \
	-param caption "odie element index" \
	-html > $@

$(DOT)/rk4tau-pages : TD := $(TARGET_DIRECTORY)
$(DOT)/rk4tau-pages : $(SCHEMA)/flat/rk4tau.rng
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc.xsl \
	-param target-directory $(TD) \
	-param toc-doc-url rk4tau-toc.html \
	-param flat-doc-url $(DOC)/overflow-html \
	-param ndx-doc-url rk4tau-index.html \
	-html
	touch $@

$(DOT)/rk4tau-toc.html : $(SCHEMA)/flat/rk4tau.rng
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc-toc.xsl \
	-param caption "rk4tau table of contents" \
	-html > $@

$(DOT)/rk4tau-index.html : $(SCHEMA)/flat/rk4tau.rng
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc-ndx.xsl \
	-param caption "rk4tau element index" \
	-html > $@

CLEAN_LIST := $(CLEAN_LIST) \
	$(DOT)/*.html \
	$(DOT)/moleculizer-pages \
	$(DOT)/moleculizer-state-pages \
	$(DOT)/odie-pages \
	$(DOT)/rk4tau-pages

PREEN_LIST := $(PREEN_LIST) $(DOT)/*~ $(DOT)/*.bak

DOT := $(call dotdot,$(DOT))
