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

DR_CURVE_PLOT := $(TDDR_DIR)/$(DMP_FILE_STEM)--$(COLUMN_HEADER)--dr.png
$(DR_CURVE_PLOT) : $(GNUPLOT_DR_CURVE_SCRIPT)
	gnuplot $< > $@

# Really, this depends on $(TDDR_HISTO_DIR)/all-stats, but that target is
# only explained in histos.mk.  Also note use of "histo-target" above, which
# is the catch-all target defined by histos.mk.
TRANSINFO_FILE := $(TDDR_DIR)/$(DMP_FILE_STEM)--$(COLUMN_HEADER)--transinfo
$(TRANSINFO_FILE) : $(TDDR_HISTO_DIR)/all-stats
	nu -load $(SCRIPT)/transinfo.scm \
	-stats-file $< \
	-transinfo-file $@

REPORT_FILE := $(TDDR_DIR)/$(DMP_FILE_STEM)--$(COLUMN_HEADER)--report
$(REPORT_FILE) : FILE_STEM := $(DMP_FILE_STEM)--$(COLUMN_HEADER)
$(REPORT_FILE) : $(TRANSINFO_FILE)
	./$(SCRIPT)/write-transinfo-report.pl $(FILE_STEM) $(VARIATIONS_FILE) $(OUTPUT_DIR)

include $(TMP)/summary-histos.mk

variations-target : $(HISTO_PLOT) \
	$(DR_CURVE_PLOT) \
	$(TRANSINFO_FILE) \
	$(REPORT_FILE)
