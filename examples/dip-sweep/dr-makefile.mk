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
	$(DD)/new-tddr.pl odie-in.xml odie-doses
	$(DD)/gen-odie-scripts.pl odie-in.xml odie-doses
	chmod a+x immed_script.sh
	chmod a+x batch_script.sh

sims-now : DD := $(DEMO_DIR)
sims-now :
	./immed_script.sh
	for dmp_file in `cd dose_0 && echo *.dmp` ; do \
	$(DD)/new-report.pl $$dmp_file odie-doses;\
	done
	cd tddrs ;\
	for tddr in *.tddr ; do \
	plt -n -y 'Concentration in moles/liter' $$tddr > $${tddr%.tddr}.png ;\
	done

# Run the dose/response simulations in the batch queue.
sims-batch :
	./batch_script.sh
