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

# For the moment, this is just testing for the new script, new-tddr.pl.

# Note that the names of these output directories have to be coordinated
# with the directory names given in odie-doses and moleculzer-doses.
OUTPUT_DIR := new-tddr.out

$(OUTPUT_DIR) :
	mkdir $@

HISTOS_DIR := $(OUTPUT_DIR)/histos

$(HISTOS_DIR) : $(OUTPUT_DIR)
	mkdir $@

# Generate randomizations of moleculizer species populations.
moleculizer-variations : moleculizer-pops
	nu -load randomize-pops.scm \
	-input $? \
	-output $@ \
	-count 10 \
	-fraction 0.33

# Set up simulations of the randomizations.  The run for each randomization
# will be a dose/response sweep, run by a makefile and a doses file.
immed_script.sh batch_script.sh : | $(OUTPUT_DIR) moleculizer-variations
	./new-tddr.pl moleculizer-in.xml moleculizer-variations
	./gen-dr-scripts.pl dr-makefile.mk \
	moleculizer-doses moleculizer-variations
	chmod +x immed_script.sh
	chmod +x batch_script.sh
	chmod +x setup_script.sh
	./setup_script.sh

sims-now : immed_script.sh
	./immed_script.sh

report : HD := $(HISTOS_DIR)
report : $(HISTOS_DIR)
	./merge-cells.pl moleculizer-variations moleculizer-doses
	./write-gp-doses-script.pl moleculizer-doses $(HD)/plot-histograms.gp $(HD)
	for DATA_FILE in $(HD)/*.data ; do \
	realbin 5 < $$DATA_FILE > $${DATA_FILE%.data}.histo ; \
	./stats.pl $${DATA_FILE%.data} < $$DATA_FILE >> $(HD)/all-stats ; \
	done
	nu -load ./transinfo.scm -stats-file $(HD)/all-stats -transinfo-file $(HD)/transinfo
	./do-gamma-densities.pl $(HD)/all-stats .
	gnuplot $(HD)/plot-histograms.gp > histograms.png
	./make-transinfo-table.pl moleculizer-variations transinfo-table.out


# Test sweep script with better values input and flexible with
# respect to target program.
# setup-odie-doses : OD := $(ODIE_OUTPUT)
# setup-odie-doses : | $(ODIE_OUTPUT)
# 	./setup-doses.pl odie-in.xml odie-doses odie
# 	for dose_dir in $(OD)/dose_* ; do \
# 	cp odie-concentrations $$dose_dir ;\
# 	nu -load randomize-rates.scm \
# 	-input $$dose_dir/odie-concentrations \
# 	-output $$dose_dir/odie-variations \
# 	-count 10 \
# 	-fraction 0.33 ;\
# 	done

clean :
	rm -rf *.out
	rm -rf batch_script.sh immed_script.sh setup_script.sh
	rm -f moleculizer-variations
	rm -f histograms.png

preen : clean
	rm -rf *.bak
	rm -rf *~
