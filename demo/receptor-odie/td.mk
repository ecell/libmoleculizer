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

DOT := $(DOT)/receptor-odie

$(DOT)/target : $(DOT)/receptor-odie.out

$(DOT)/receptor-odie.out : $(DOT)/receptor-odie.xml
	mkdir $@
	cp $? $@
	cd $@ \
	&&  odie < receptor-odie.xml \
	&& plot-dmp-files-as-cnc

# This automatically generated odie file has to be edited to get the correct
# stop time, epsilons, etc.
$(DOT)/receptor-odie.xml.out : $(DOT)/receptor-old.out
	state2odie $?/param-dump.xml $@

# Generate a dump file from a version of the old receptor demo. The old
# receptor demo is altered so that alpha is put in at the beginning, instead
# of in the middle of the simulation.
$(DOT)/receptor-old.out : $(DOT)/receptor-old.xml
	mkdir $@
	cp $? $@
	cd $@ \
	&& moleculizer < receptor-old.xml \
	&& parametrizer stateDump < receptor-old.xml > param-dump.xml \
	&& plot-dmp-files

# For redoing individual demos several times.
$(DOT)/clean : DT := $(DOT)
$(DOT)/clean :
	rm -rf $(DT)/*.out

CLEAN_LIST := $(CLEAN_LIST) $(DOT)/*.out

PREEN_LIST := $(PREEN_LIST) $(DOT)/*~ $(DOT)/*.bak

DOT := $(call dotdot,$(DOT))
