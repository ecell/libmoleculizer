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

# This makefile controls a single dose/response simulation, the second
# tier of this two-tiered sweep.

DEMO_DIR := ../..

# Set up the subdirectories and scripts for the dose/response simulations.
setup : DD := $(DEMO_DIR)
setup :
	$(DD)/new-tddr.pl moleculizer-in.xml moleculizer-doses
	$(DD)/gen-mzr-scripts.pl moleculizer-in.xml moleculizer-doses
	chmod a+x immed_script.sh
	chmod a+x batch_script.sh

# Run the dose/response simulations now.  For testing, adding on all the
# commands from report, too.
#
# I note that generating the gnuplot script here is entirely redundant.
# This script is the same everywhere and it should be written once and
# for all at the top level.
sims-now : DD := $(DEMO_DIR)
sims-now :
	./immed_script.sh
	for dmp_file in `cd dose_0 && echo *.dmp` ; do \
	$(DD)/new-report.pl $$dmp_file moleculizer-doses;\
	done
	cd tddrs ;\
	for tddr in *.tddr ; do \
	plt -n $$tddr > $${tddr%.tddr}.png ;\
	done
	mkdir histos
	$(DD)/tail-histos.pl tddrs/active-kinases--active-kin2.tddr \
	200 5 histos
	$(DD)/write-gp-doses-script.pl moleculizer-doses histos/plot-histograms.gp histos
	for DATA_FILE in histos/*.data ; do \
	realbin 5 < $$DATA_FILE > $${DATA_FILE%.data}.histo ; \
	$(DD)/stats.pl $${DATA_FILE%.data} < $$DATA_FILE >> histos/all-stats ; \
	done
	nu -load $(DD)/transinfo.scm -stats-file histos/all-stats -transinfo-file histos/transinfo
	$(DD)/do-gamma-densities.pl histos/all-stats $(DD)
	gnuplot histos/plot-histograms.gp > histograms.png

# Run the dose/response simulations in the batch queue.
sims-batch :
	./batch_script.sh

# Generate tddr report files and do simple plots.
report : DD := $(DEMO_DIR)
report :
	for dmp_file in `cd sim_0 && echo *.dmp` ; do \
	$(DD)/new-report.pl $$dmp_file moleculizer-doses;\
	done
	cd tddrs ;\
	for tddr in *.tddr ; do \
	plt -n $$tddr > $${tddr%.tddr}.png ;\
	done
	mkdir histos
	$(DD)/tail-histos.pl tddrs/active-kinases--active-kin2.tddr \
	200 5 histos
	nu -load transinfo.scm -stats-file histos/all-stats
	$(DD)/do-gamma-densities.sh
	gnuplot plot-histograms.gp > histograms.png

