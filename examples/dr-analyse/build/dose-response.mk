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

$(DR_OUTPUT_DIR) : $(OUTPUT_DIR)
	mkdir $@

dose-response-dirs : | $(DR_OUTPUT_DIR)

TDDR_DIR := $(DR_OUTPUT_DIR)/tddrs

$(TDDR_DIR) : | $(DR_OUTPUT_DIR)
	mkdir $@

TEMPLATE_INPUT := $(DR_OUTPUT_DIR)/moleculizer-in.xml
DOSES_FILE := ./moleculizer-doses

$(DR_OUTPUT_DIR)/doses-setup : DSES_FILE := $(DOSES_FILE)
$(DR_OUTPUT_DIR)/doses-setup : TPLT_INPUT := $(TEMPLATE_INPUT)
$(DR_OUTPUT_DIR)/doses-setup : variations-setup $(DOSES_FILE) \
	$(DR_OUTPUT_DIR)/dose_dirs
	$(SCRIPT)/do-substitutions.pl $(TPLT_INPUT) \
	$(DSES_FILE) $(@D) moleculizer-in.xml

# Probably will also need a fig version of this plot.
GNUPLOT_DORA_PLOT_SCRIPT := $(DR_OUTPUT_DIR)/dora-plot-script.gp
$(GNUPLOT_DORA_PLOT_SCRIPT) : TDDRD := $(TDDR_DIR)
$(GNUPLOT_DORA_PLOT_SCRIPT) : | $(DR_OUTPUT_DIR)
	$(SCRIPT)/write-gp-dora-plot-script.pl \
	$(TDDRD) $(RESPONSES_FILE) $@

DORA_PNG_PLOT := $(DR_OUTPUT_DIR)/dora-plot.png
$(DORA_PNG_PLOT) : $(GNUPLOT_DORA_PLOT_SCRIPT) $(DR_OUTPUT_DIR)/tddr-target
	gnuplot $< > $@

# N is for "normalized."
GNUPLOT_NDORA_PLOT_SCRIPT := $(DR_OUTPUT_DIR)/ndora-plot-script.gp
$(GNUPLOT_NDORA_PLOT_SCRIPT) : TDDRD := $(TDDR_DIR)
$(GNUPLOT_NDORA_PLOT_SCRIPT) : | $(DR_OUTPUT_DIR)
	$(SCRIPT)/write-gp-ndora-plot-script.pl \
	$(TDDRD) $(RESPONSES_FILE) $@

NDORA_PNG_PLOT := $(DR_OUTPUT_DIR)/normalized-dora-plot.png
$(NDORA_PNG_PLOT) : $(GNUPLOT_NDORA_PLOT_SCRIPT) $(DR_OUTPUT_DIR)/tddr-target
	gnuplot $< > $@

HILL_REPORT := $(DR_OUTPUT_DIR)/hill-report.txt
$(HILL_REPORT) : TDDRD := $(TDDR_DIR)
$(HILL_REPORT) : | $(DR_OUTPUT_DIR)
	$(SCRIPT)/write-hill-report.pl \
	$(TDDRD) $(RESPONSES_FILE) $@

variations-target : $(DORA_PNG_PLOT) \
	$(NDORA_PNG_PLOT) \
	$(HILL_REPORT)

-include $(TMP)/doses.mk
-include $(TMP)/responses.mk
