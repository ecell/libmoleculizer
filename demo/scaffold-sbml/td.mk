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

DOT := $(DOT)/scaffold-sbml

$(DOT)/target : $(DOT)/sbml-pretty.txt $(DOT)/scaffold-state-pretty.txt

# Pretty-prints sbml output.
$(DOT)/sbml-pretty.txt : $(DOT)/sbml.out
	xmlpretty $? $@

# Pretty-prints state dump output.
$(DOT)/scaffold-state-pretty.txt : $(DOT)/stateDump
	xmlpretty $? $@

# Converts state dump to sbml.
$(DOT)/sbml.out : $(DOT)/stateDump
	state2sbml $? $@

# Uncompresses state dump from alpha pathway simulation done off-line.
$(DOT)/stateDump : DT := $(DOT)
$(DOT)/stateDump : $(DOT)/stateDump.tgz
	cd $(DT) && tar -xvzf stateDump.tgz

# Runs alpha pathway simulation to generate stateDump, which
# it then compresses.
$(DOT)/stateDump.tgz : DT := $(DOT)
$(DOT)/stateDump.tgz :
	mkdir $(DT)/scaffold.out
	cp $(DT)/scaffold.xml $(DT)/scaffold.out
	cd $(DT)/scaffold.out \
	&& moleculizer < scaffold.xml \
	&& plot-dmp-files \
	&& tar -cvzf stateDump.tgz stateDump
	cp $(DT)/scaffold.out/stateDump.tgz $(DT)

# For redoing individual demos several times.
$(DOT)/clean : DT := $(DOT)
$(DOT)/clean :
	rm -rf $(DT)/*.out \
	$(DT)/stateDump \
	$(DT)/sbml-pretty.txt \
	$(DT)/scaffold-state-pretty.txt

CLEAN_LIST := $(CLEAN_LIST) \
	$(DOT)/*.out \
	$(DOT)/stateDump \
	$(DOT)/sbml-pretty.txt \
	$(DOT)/scaffold-state-pretty.txt

PREEN_LIST := $(PREEN_LIST) $(DOT)/*~ $(DOT)/*.bak

DOT := $(call dotdot,$(DOT))
