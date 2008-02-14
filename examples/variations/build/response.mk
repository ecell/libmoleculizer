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

TDDR_FILE := $(TDDR_DIR)/$(DMP_FILE_STEM)--$(COLUMN_HEADER).tddr

# cut-dmp-column.pl leaves useful intermediate output in the TDDR_HISTO_DIR,
# so that this TDDR_FILE is a "prerequisite" for things derived from the
# intermediate output.
$(TDDR_FILE) : DROD := $(DR_OUTPUT_DIR)
$(TDDR_FILE) : STEM := $(DMP_FILE_STEM)
$(TDDR_FILE) : HDR := $(COLUMN_HEADER)
$(TDDR_FILE) : TDDRD := $(TDDR_DIR)
$(TDDR_FILE) :  doses | $(TDDR_HISTO_DIR) $(TDDR_DIR)
	./$(SCRIPT)/cut-dmp-column.pl $(STEM) $(HDR) \
	$(DOSES_FILE) $(DROD) $(TDDRD)

TDDR_PNG := $(TDDR_DIR)/$(DMP_FILE_STEM)--$(COLUMN_HEADER).png

$(TDDR_PNG) : $(TDDR_FILE)
	plt -n $< > $@

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

TRANSINFO_FILE := $(TDDR_DIR)/$(DMP_FILE_STEM)--$(COLUMN_HEADER)--transinfo
$(TRANSINFO_FILE) : $(TDDR_HISTO_DIR)/all-stats
	nu -load $(SCRIPT)/transinfo.scm \
	-stats-file $< \
	-transinfo-file $@

-include $(TMP)/histos.mk

$(DR_OUTPUT_DIR)/tddr-target : $(TDDR_PNG) \
	$(HISTO_PLOT) \
	$(DR_CURVE_PLOT) \
	$(TRANSINFO_FILE)
