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
$(TDDR_FILE) : DROD := $(OUTPUT_DIR)
$(TDDR_FILE) : STEM := $(DMP_FILE_STEM)
$(TDDR_FILE) : HDR := $(COLUMN_HEADER)
$(TDDR_FILE) : TDDRD := $(TDDR_DIR)
$(TDDR_FILE) :  $(OUTPUT_DIR)/doses | $(TDDR_HISTO_DIR) $(TDDR_DIR)
	./$(SCRIPT)/cut-dmp-column.pl $(STEM) $(HDR) \
	$(DOSES_FILE) $(DROD) $(TDDRD)

# TRANSINFO_FILE := $(TDDR_DIR)/$(DMP_FILE_STEM)--$(COLUMN_HEADER)--transinfo
# $(TRANSINFO_FILE) : $(TDDR_HISTO_DIR)/all-stats
# 	nu -load $(SCRIPT)/transinfo.scm \
# 	-stats-file $< \
# 	-transinfo-file $@

# Use new infocap executable to calculate information capacity and tabulate
# optimal dose density function.

# Probably will want to use one of the other transinformation-computing
# commands.  Presently, this is using the derivative of the dose/response
# curve as dose density, so derivTI should be more or less equivalent to this.
# Also have gammaTI and uniformTI available.

# Formerly used this command, which does optimization.  
# infocap -ms 1.0e3 -sf $< -of $@ -df $(DF)
TRANSINFO_FILE := $(TDDR_DIR)/$(DMP_FILE_STEM)--$(COLUMN_HEADER)--transinfo
$(TRANSINFO_FILE) : DF := $(TDDR_HISTO_DIR)/optimal-density
$(TRANSINFO_FILE) : $(TDDR_HISTO_DIR)/all-stats
	infocap -sf $< -dt $@ -df $(DF)

-include $(TMP)/histos.mk

$(OUTPUT_DIR)/all-histos : $(TDDR_HISTO_DIR)/histo-target

target : $(TRANSINFO_FILE)
