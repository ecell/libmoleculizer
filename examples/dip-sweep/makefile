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

# Note that the names of these output directories have to be coordinated
# with the directory names given in odie-doses and moleculzer-doses.
OUTPUT_DIR := sweep.out

$(OUTPUT_DIR) :
	mkdir $@

# Set up simulations of the randomizations.  The run for each randomization
# will be a dose/response sweep, run by a makefile and a doses file.
immed_script.sh batch_script.sh : | $(OUTPUT_DIR) feedback-rates
	./new-tddr.pl odie-in.xml feedback-rates
	./gen-dr-scripts.pl dr-makefile.mk \
	odie-doses feedback-rates
	chmod +x immed_script.sh
	chmod +x batch_script.sh
	chmod +x setup_script.sh
	./setup_script.sh

sims-now : immed_script.sh
	./immed_script.sh

clean :
	rm -rf *.out
	rm -f batch_script.sh immed_script.sh setup_script.sh

preen : clean
	rm -rf *.bak
	rm -rf *~
