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

DOT := $(DOT)/demo

SHORT_DEMOS := 	bndKinase \
	continuator \
	heinrich \
	heinrich-odie \
	heinrich-rk4tau \
	heinrich-sbml \
	heinrich-state \
	kinase \
	kinase-extrap \
	kinase-odie \
	kinase-sbml \
	modify \
	omniKinase \
	query-allostery \
	receptor \
	receptor-nx-extrap \
	receptor-nx-sst2 \
	receptor-odie \
	receptor-sbml \
	receptor-state \
	receptor-structure \
	scaffold-odie \
	scaffold-sbml \
	scaffold-state \
	simple \
	simple-extrap \
	simple-odie \
	simple-rk4tau \
	simple-sbml \
	simple-stoch \
	small-mol

LONG_DEMOS := scaffold \
	scaffold-basal \
	scaffold-extrap

DEMOS := $(SHORT_DEMOS) $(LONG_DEMOS)

$(DOT)/short : $(addprefix $(DOT)/,$(addsuffix /target,$(SHORT_DEMOS)))

$(DOT)/long : $(addprefix $(DOT)/,$(addsuffix /target,$(LONG_DEMOS)))

$(DOT)/all : $(addprefix $(DOT)/,$(addsuffix /target,$(DEMOS)))

$(DOT)/clean : $(addprefix $(DOT)/,$(addsuffix /clean,$(DEMOS)))

PREEN_LIST := $(PREEN_LIST) $(DOT)/*~

include $(addprefix $(DOT)/,$(addsuffix /td.mk,$(DEMOS)))

DOT := $(call dotdot,$(DOT))

