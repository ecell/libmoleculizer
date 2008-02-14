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

$(DR_OUTPUT_DIR) : $(VAR_SUITE_OUTPUT_DIR)
	mkdir $@

# This target exists so that all the directories can be created before the
# perl script that generates the substituted moleculizer input files runs.
$(VAR_SUITE_OUTPUT_DIR)/dose-response-dirs : | $(DR_OUTPUT_DIR)

TDDR_DIR := $(DR_OUTPUT_DIR)/tddrs

$(TDDR_DIR) : | $(DR_OUTPUT_DIR)
	mkdir $@

TEMPLATE_INPUT := $(DR_OUTPUT_DIR)/moleculizer-in.xml
DOSES_FILE := ./moleculizer-doses

# Setup of the doses input scripts.
# 
# The target $(VAR_SUITE_OUTPUT_DIR)/variations-setup ensures that the
# directories where the input scripts go exist before they are put there.
$(DR_OUTPUT_DIR)/doses-setup : DSES_FILE := $(DOSES_FILE)
$(DR_OUTPUT_DIR)/doses-setup : TPLT_INPUT := $(TEMPLATE_INPUT)
$(DR_OUTPUT_DIR)/doses-setup : RT_DIR := $(DR_OUTPUT_DIR)
$(DR_OUTPUT_DIR)/doses-setup : $(VAR_SUITE_OUTPUT_DIR)/variations-setup \
	$(TEMPLATE_INPUT) \
	$(DOSES_FILE) \
	$(DR_OUTPUT_DIR)/dose_dirs
	$(SCRIPT)/do-substitutions.pl $(TPLT_INPUT) \
	$(DSES_FILE) $(RT_DIR) moleculizer-in.xml

# The summary dose/response curve plot, where all the dose/response
# curves are plotted on one plot.  This is for eyeballing the "xerox effect".
# 
# The $(DR_OUTPUT_DIR)/all-histos target guarantees that the stats files from
# which this plot is created all exist.
$(DR_OUTPUT_DIR)/dr-curve-summary.gp : TDDRD := $(TDDR_DIR)
$(DR_OUTPUT_DIR)/dr-curve-summary.gp : $(DR_OUTPUT_DIR)/all-histos
	$(SCRIPT)/write-summary-dr-curve-script.pl \
	$(TDDRD) $(RESPONSES_FILE) $@

$(DR_OUTPUT_DIR)/dr-curve-summary.png : DOD := $(DR_OUTPUT_DIR)
$(DR_OUTPUT_DIR)/dr-curve-summary.png : $(DR_OUTPUT_DIR)/dr-curve-summary.gp
	gnuplot $(DOD)/dr-curve-summary.gp > $@

$(DR_OUTPUT_DIR)/fig-dr-curve-summary.gp : TDDRD := $(TDDR_DIR)
$(DR_OUTPUT_DIR)/fig-dr-curve-summary.gp : $(DR_OUTPUT_DIR)/all-histos
	$(SCRIPT)/write-fig-summary-dr-curve-script.pl \
	$(TDDRD) $(RESPONSES_FILE) $@

$(DR_OUTPUT_DIR)/dr-curve-summary.fig : DOD := $(DR_OUTPUT_DIR)
$(DR_OUTPUT_DIR)/dr-curve-summary.fig : $(DR_OUTPUT_DIR)/fig-dr-curve-summary.gp
	gnuplot $(DOD)/fig-dr-curve-summary.gp > $@

$(DR_OUTPUT_DIR)/dora-score : TDDRD := $(TDDR_DIR)
$(DR_OUTPUT_DIR)/dora-score : $(DR_OUTPUT_DIR)/all-histos
	$(SCRIPT)/dora-score.pl $(TDDRD) $(RESPONSES_FILE)

$(VAR_SUITE_OUTPUT_DIR)/variations : \
	$(DR_OUTPUT_DIR)/dr-curve-summary.png \
	$(DR_OUTPUT_DIR)/dr-curve-summary.fig \
	$(DR_OUTPUT_DIR)/dora-score

-include $(TMP)/doses.mk
-include $(TMP)/responses.mk
