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

TDDR_HISTO_DIR := $(TDDR_DIR)/$(DMP_FILE_STEM)--$(COLUMN_HEADER).histos
$(TDDR_HISTO_DIR) : | $(TDDR_DIR)
	mkdir $@

GNUPLOT_HISTO_SCRIPT := $(TDDR_HISTO_DIR)/histos.gp
$(GNUPLOT_HISTO_SCRIPT) : THD := $(TDDR_HISTO_DIR)
$(GNUPLOT_HISTO_SCRIPT) : | $(TDDR_HISTO_DIR)
	./$(SCRIPT)/write-gp-histos-script.pl $(DOSES_FILE) $@ $(THD)

HISTO_PLOT := $(TDDR_DIR)/$(DMP_FILE_STEM)--$(COLUMN_HEADER)--histo.png
$(HISTO_PLOT) : $(GNUPLOT_HISTO_SCRIPT) $(TDDR_HISTO_DIR)/histo-target
	gnuplot $< > $@

GNUPLOT_DR_CURVE_SCRIPT := $(TDDR_HISTO_DIR)/dr_curve.gp
$(GNUPLOT_DR_CURVE_SCRIPT) : THD := $(TDDR_HISTO_DIR)
$(GNUPLOT_DR_CURVE_SCRIPT) : | $(TDDR_HISTO_DIR)
	./$(SCRIPT)/write-gp-dr-curve-script.pl $(THD) $@

GNUPLOT_FIG_DR_CURVE_SCRIPT := $(TDDR_HISTO_DIR)/fig_dr_curve.gp
$(GNUPLOT_FIG_DR_CURVE_SCRIPT) : THD := $(TDDR_HISTO_DIR)
$(GNUPLOT_FIG_DR_CURVE_SCRIPT) : | $(TDDR_HISTO_DIR)
	./$(SCRIPT)/write-fig-dr-curve-script.pl $(THD) $@

DR_CURVE_PLOT := $(TDDR_DIR)/$(DMP_FILE_STEM)--$(COLUMN_HEADER)--dr.png
$(DR_CURVE_PLOT) : $(GNUPLOT_DR_CURVE_SCRIPT)
	gnuplot $< > $@

DR_CURVE_FIG := $(TDDR_DIR)/$(DMP_FILE_STEM)--$(COLUMN_HEADER)--dr.fig
$(DR_CURVE_FIG) : $(GNUPLOT_FIG_DR_CURVE_SCRIPT)
	gnuplot $< > $@

# Really, this depends on $(TDDR_HISTO_DIR)/all-stats, but that target is
# only explained in histos.mk.  Also note use of "histo-target" above, which
# is the catch-all target defined by histos.mk.

# TRANSINFO_FILE := $(TDDR_DIR)/$(DMP_FILE_STEM)--$(COLUMN_HEADER)--transinfo
# $(TRANSINFO_FILE) : $(TDDR_HISTO_DIR)/all-stats
# 	nu -load $(SCRIPT)/transinfo.scm \
# 	-stats-file $< \
# 	-transinfo-file $@

TRANSINFO_FILE := $(TDDR_DIR)/$(DMP_FILE_STEM)--$(COLUMN_HEADER)--transinfo
$(TRANSINFO_FILE) : $(TDDR_HISTO_DIR)/all-stats
	infocap -sf $< -of $@

REPORT_FILE := $(TDDR_DIR)/$(DMP_FILE_STEM)--$(COLUMN_HEADER)--report
$(REPORT_FILE) : FILE_STEM := $(DMP_FILE_STEM)--$(COLUMN_HEADER)
$(REPORT_FILE) : VSOD := $(VAR_SUITE_OUTPUT_DIR)
$(REPORT_FILE) : $(TRANSINFO_FILE) $(VAR_SUITE_OUTPUT_DIR)/variations
	./$(SCRIPT)/write-transinfo-report.pl $(FILE_STEM) $(VARIATIONS_FILE) $(VSOD)

SCATTERPLOT_FILE := $(TDDR_DIR)/$(DMP_FILE_STEM)--$(COLUMN_HEADER)--scatterplot
$(SCATTERPLOT_FILE) : FILE_STEM := $(DMP_FILE_STEM)--$(COLUMN_HEADER)
$(SCATTERPLOT_FILE) : VSOD := $(VAR_SUITE_OUTPUT_DIR)
$(SCATTERPLOT_FILE) : $(TRANSINFO_FILE) $(VAR_SUITE_OUTPUT_DIR)/variations
	./$(SCRIPT)/write-scatterplot.pl $(FILE_STEM) $(VARIATIONS_FILE) $(VSOD)

GNUPLOT_SCATTERPLOT_SCRIPT :=  $(TDDR_DIR)/$(DMP_FILE_STEM)--$(COLUMN_HEADER)--scatterplot.gp
$(GNUPLOT_SCATTERPLOT_SCRIPT) : SLT := $(SCATTERPLOT_FILE)
$(GNUPLOT_SCATTERPLOT_SCRIPT) : | $(TDDR_DIR)
	./$(SCRIPT)/write-gp-scatterplot-script.pl $(SLT) $@

SCATTERPLOT_PLOT := $(TDDR_DIR)/$(DMP_FILE_STEM)--$(COLUMN_HEADER)--scatterplot.png
$(SCATTERPLOT_PLOT) : $(GNUPLOT_SCATTERPLOT_SCRIPT) $(SCATTERPLOT_FILE)
	gnuplot $< > $@

include $(TMP)/summary-histos.mk

$(OUTPUT_DIR)/feedbacks : $(HISTO_PLOT) \
	$(DR_CURVE_PLOT) \
	$(DR_CURVE_FIG) \
	$(TRANSINFO_FILE) \
	$(REPORT_FILE) \
	$(SCATTERPLOT_PLOT)
