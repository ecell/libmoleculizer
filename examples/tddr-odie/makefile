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

# Generates and runs simulations now (not in batch queue), then
# generates tddr report.
report-now : sims-now
	$(MAKE) report

# Sets up simulations for tddr, writes batch script and immediate
# simulation script.
setup :
	tddr-float-setup-odie.pl in.xml values sim
	chmod a+x ./batch_script.sh
	chmod a+x ./immed_script.sh

# Run simulations in batch queue.
sims-batch : setup
	./batch_script.sh 2>&1 | tee batch_script.log

# Run simulations now.
sims-now : setup
	./immed_script.sh 2>&1 | tee immed_script.log

# Generate tddr report files and do simple plots.
report :
	for dmp_file in `cd sim_0 && echo *.dmp` ; do \
	tddr-report.pl $$dmp_file values sim ;\
	done
	cd tddrs ;\
	for tddr in *.tddr ; do \
	plt -n $$tddr > $${tddr%.tddr}.png ;\
	done

# General clean-up.
clean :
	rm -rf sim_*
	rm -f batch_script.sh
	rm -f batch_script.log
	rm -f immed_script.sh
	rm -f immed_script.log
	rm -rf tddrs

preen : clean
	rm -f *~ *.bak
