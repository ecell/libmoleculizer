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

DOT := $(DOT)/scripts

SCRIPTS := altSeedString \
	concat.pl \
	getSeedString \
	gpFunctions.pm \
	mzr-dispatch-task \
	plot-cnc-files \
	plot-dmp-files \
	plot-dmp-files-as-cnc \
	plt \
	ready-targets.bash \
	state2odie \
	state2rk4tau \
	state2sbml \
	vary-param \

BIN_SCRIPTS := $(addprefix $(BIN)/,$(SCRIPTS))

$(DOT)/target : $(BIN_SCRIPTS)

# Here avoiding the obvious application of 'eval' which would work fine
# here, since the strings are relatively short.
$(BIN)/altSeedString : $(DOT)/altSeedString
	cp $< $@

$(BIN)/concat.pl : $(DOT)/concat.pl
	cp $< $@

$(BIN)/getSeedString : $(DOT)/getSeedString
	cp $< $@

$(BIN)/gpFunctions.pm : $(DOT)/gpFunctions.pm
	cp $< $@

$(BIN)/mzr-dispatch-task : $(DOT)/mzr-dispatch-task
	cp $< $@

$(BIN)/plot-cnc-files : $(DOT)/plot-cnc-files
	cp $< $@

$(BIN)/plot-dmp-files : $(DOT)/plot-dmp-files
	cp $< $@

$(BIN)/plot-dmp-files-as-cnc : $(DOT)/plot-dmp-files-as-cnc
	cp $< $@

$(BIN)/plt : $(DOT)/plt
	cp $< $@

$(BIN)/ready-targets.bash : $(DOT)/ready-targets.bash
	cp $< $@

$(BIN)/state2odie : $(DOT)/state2odie
	cp $< $@

$(BIN)/state2rk4tau : $(DOT)/state2rk4tau
	cp $< $@

$(BIN)/state2sbml : $(DOT)/state2sbml
	cp $< $@

$(BIN)/vary-param : $(DOT)/vary-param
	cp $< $@

CLEAN_LIST := $(CLEAN_LIST) $(BIN_SCRIPTS)

PREEN_LIST := $(PREEN_LIST) $(DOT)/*~

DOT := $(call dotdot,$(DOT))
