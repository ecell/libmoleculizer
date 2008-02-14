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

OFLOW_FILES := altBlobGraph.jpg \
	cptzr.html \
	cptzr-flowchart.jpg \
	gpa-exchange-gen.html \
	hetero-hydrolysis-gen.html \
	install-guide.html \
	moleculizer-input.html \
	mols.html \
	more-gpa-exchange-gen.html \
	mzr-variation-actions.html \
	mzr-variation-upload.html \
	nucleotide-exchange-gen.html \
	reparam-upload.html \
	rk4tau-variation-actions.html \
	rk4tau-variation-upload.html \
	state-reparam-actions.html \
	state-reparam-upload.html \
	stat-stream-ref.html \
	user-guide.html \
	welcome.html

OFLOW_SRC_DIR := $(DOCS)/overflow-html

OFLOW_SRC_FILES := $(addprefix $(OFLOW_SRC_DIR)/,$(OFLOW_FILES))

OFLOW_TGT_DIR := $(DOC)/overflow-html

OFLOW_TGT_FILES := $(addprefix $(OFLOW_TGT_DIR)/,$(OFLOW_FILES))

$(OFLOW_TGT_DIR) : | $(DOC)
	mkdir $@

# Copy hand-generated html file from the source area to doc.
$(OFLOW_TGT_FILES) : $(OFLOW_TGT_DIR)/% : $(OFLOW_SRC_DIR)/% | $(OFLOW_TGT_DIR)
	cp $< $@

# Add overflow html to overall doc target.
usr-doc : $(OFLOW_TGT_FILES)

PREEN_LIST += $(DOCS)/overflow-html/*.bak




