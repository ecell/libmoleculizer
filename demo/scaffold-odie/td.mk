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

DOT := $(DOT)/scaffold-odie

$(DOT)/target : $(DOT)/scaffold-odie.out

# The full scaffold demo using odie.
$(DOT)/scaffold-odie.out : $(DOT)/scaffold-odie.xml
	mkdir $@
	cp $? $@
	cd $@ \
	&& odie < scaffold-odie.xml \
	&& plot-dmp-files-as-cnc

# This file has to be edited, for now. Perhaps the additional material should
# go in an XML file, or in arguments to the XSLT transformation?
$(DOT)/scaffold-odie.xml.out : $(DOT)/scaffold.out
	state2odie $?/param-dump.xml $@

# The moleculizer simulation of the alpha pathway down through
# the scaffold complex.
$(DOT)/scaffold.out : $(DOT)/scaffold.xml
	mkdir $@
	cp $? $@
	cd $@ \
	&& moleculizer < scaffold.xml \
	&& parametrizer state-dump.xml < scaffold.xml > param-dump.xml \
	&& plot-dmp-files

# For redoing individual demos several times.
$(DOT)/clean : DT := $(DOT)
$(DOT)/clean :
	rm -rf $(DT)/*.out

CLEAN_LIST := $(CLEAN_LIST) $(DOT)/*.out

PREEN_LIST := $(PREEN_LIST) $(DOT)/*~ $(DOT)/*.bak

DOT := $(call dotdot,$(DOT))
