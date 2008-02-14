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

STATIC_DOC_DIR := $(DOC)/static-doc

$(STATIC_DOC_DIR) : | $(DOC)
	mkdir $@

# This absolute file path is needed to handle a special property of Xalan
# XSLT processing: paths used in the redirect feature are normally taken
# to be relative the file of the input file.  In these cases, that would
# be somewhere off in XML land.
#
# This absolute path doesn't actually go into any of the pages themselves,
# which are linked to each other by relative links.
# 
# CURDIR is defined automatically by make and works well with $(MAKE) -C.
DOC_TREE_DIR := $(CURDIR)/$(STATIC_DOC_DIR)

# Here, it looks like flat-doc-url has to be the relative path from the main
# bolus of pages to the overflow pages.  Previously, I had
# $(DOC)/overflow-html, which doesn't work at all.
$(STATIC_DOC_DIR)/.moleculizer-pages-target : TD := $(DOC_TREE_DIR)
$(STATIC_DOC_DIR)/.moleculizer-pages-target : $(FLATTENED_DIR)/moleculizer.rng | $(STATIC_DOC_DIR)
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc.xsl \
	-param target-directory $(TD) \
	-param toc-doc-url moleculizer-toc.html \
	-param flat-doc-url '../overflow-html' \
	-param ndx-doc-url moleculizer-index.html \
	-html
	touch $@

$(STATIC_DOC_DIR)/moleculizer-toc.html : $(FLATTENED_DIR)/moleculizer.rng | $(STATIC_DOC_DIR)
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc-toc.xsl \
	-param caption "moleculizer-input table of contents" \
	-html > $@

$(STATIC_DOC_DIR)/moleculizer-index.html : $(FLATTENED_DIR)/moleculizer.rng | $(STATIC_DOC_DIR)
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc-ndx.xsl \
	-param caption "moleculizer-input element index" \
	-html > $@

###################

# Consider regularizing the name of the schema here.  Could the rules for the
# various schema be combined into one pattern-based rule?

$(STATIC_DOC_DIR)/.cpt-pages-target : TD := $(DOC_TREE_DIR)
$(STATIC_DOC_DIR)/.cpt-pages-target : $(FLATTENED_DIR)/compartment.rng | $(STATIC_DOC_DIR)
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc.xsl \
	-param target-directory $(TD) \
	-param toc-doc-url cpt-toc.html \
	-param flat-doc-url '../overflow-html' \
	-param ndx-doc-url cpt-index.html \
	-html
	touch $@

$(STATIC_DOC_DIR)/cpt-toc.html : $(FLATTENED_DIR)/compartment.rng | $(STATIC_DOC_DIR)
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc-toc.xsl \
	-param caption "compartment-model table of contents" \
	-html > $@

$(STATIC_DOC_DIR)/cpt-index.html : $(FLATTENED_DIR)/compartment.rng | $(STATIC_DOC_DIR)
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc-ndx.xsl \
	-param caption "compartment-model element index" \
	-html > $@

###################

$(STATIC_DOC_DIR)/.moleculizer-state-pages-target : TD := $(DOC_TREE_DIR)
$(STATIC_DOC_DIR)/.moleculizer-state-pages-target : $(FLATTENED_DIR)/mzrState.rng | $(STATIC_DOC_DIR)
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc.xsl \
	-param target-directory $(TD) \
	-param toc-doc-url moleculizer-state-toc.html \
	-param flat-doc-url '../overflow-html' \
	-param ndx-doc-url moleculizer-state-index.html \
	-html
	touch $@

$(STATIC_DOC_DIR)/moleculizer-state-toc.html : $(FLATTENED_DIR)/mzrState.rng | $(STATIC_DOC_DIR)
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc-toc.xsl \
	-param caption "moleculizer-state table of contents" \
	-html > $@

$(STATIC_DOC_DIR)/moleculizer-state-index.html : $(FLATTENED_DIR)/mzrState.rng | $(STATIC_DOC_DIR)
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc-ndx.xsl \
	-param caption "moleculizer-input element index" \
	-html > $@

###################

$(STATIC_DOC_DIR)/.odie-pages-target : TD := $(DOC_TREE_DIR)
$(STATIC_DOC_DIR)/.odie-pages-target : $(FLATTENED_DIR)/odie.rng | $(STATIC_DOC_DIR)
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc.xsl \
	-param target-directory $(TD) \
	-param toc-doc-url odie-toc.html \
	-param flat-doc-url '../overflow-html' \
	-param ndx-doc-url odie-index.html \
	-html
	touch $@

$(STATIC_DOC_DIR)/odie-toc.html : $(FLATTENED_DIR)/odie.rng | $(STATIC_DOC_DIR)
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc-toc.xsl \
	-param caption "odie table of contents" \
	-html > $@

$(STATIC_DOC_DIR)/odie-index.html : $(FLATTENED_DIR)/odie.rng | $(STATIC_DOC_DIR)
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc-ndx.xsl \
	-param caption "odie element index" \
	-html > $@

###################

$(STATIC_DOC_DIR)/.rk4tau-pages-target : TD := $(DOC_TREE_DIR)
$(STATIC_DOC_DIR)/.rk4tau-pages-target : $(FLATTENED_DIR)/rk4tau.rng | $(STATIC_DOC_DIR)
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc.xsl \
	-param target-directory $(TD) \
	-param toc-doc-url rk4tau-toc.html \
	-param flat-doc-url '../overflow-html' \
	-param ndx-doc-url rk4tau-index.html \
	-html
	touch $@

$(STATIC_DOC_DIR)/rk4tau-toc.html : $(FLATTENED_DIR)/rk4tau.rng | $(STATIC_DOC_DIR)
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc-toc.xsl \
	-param caption "rk4tau table of contents" \
	-html > $@

$(STATIC_DOC_DIR)/rk4tau-index.html : $(FLATTENED_DIR)/rk4tau.rng | $(STATIC_DOC_DIR)
	java org.apache.xalan.xslt.Process \
	-in $< \
	-xsl $(XSL)/static-doc-ndx.xsl \
	-param caption "rk4tau element index" \
	-html > $@

###################

# Add the files and directories full of files to the general doc target.
usr-doc : $(STATIC_DOC_DIR)/.moleculizer-pages-target \
	$(STATIC_DOC_DIR)/moleculizer-toc.html \
	$(STATIC_DOC_DIR)/moleculizer-index.html \
	$(STATIC_DOC_DIR)/.cpt-pages-target \
	$(STATIC_DOC_DIR)/cpt-toc.html \
	$(STATIC_DOC_DIR)/cpt-index.html \
	$(STATIC_DOC_DIR)/.moleculizer-state-pages-target \
	$(STATIC_DOC_DIR)/moleculizer-state-toc.html \
	$(STATIC_DOC_DIR)/moleculizer-state-index.html \
	$(STATIC_DOC_DIR)/.odie-pages-target \
	$(STATIC_DOC_DIR)/odie-toc.html \
	$(STATIC_DOC_DIR)/odie-index.html \
	$(STATIC_DOC_DIR)/.rk4tau-pages-target \
	$(STATIC_DOC_DIR)/rk4tau-toc.html \
	$(STATIC_DOC_DIR)/rk4tau-index.html
