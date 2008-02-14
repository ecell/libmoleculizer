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

DEMO_NAMES := continuator \
	cpt-omniKinase \
	cpt-receptor \
	cpt-receptor10 \
	cpt-receptor-sphere \
	cpt-simple \
	cpt-simple-stoch \
	cpt-uniMolPtase \
	feedback \
	feedback-odie \
	heinrich \
	heinrich-odie \
	heinrich-rk4tau \
	heinrich-sbml \
	heinrich-state \
	omniKinase \
	omniPtase \
	omniReceptor \
	query-allostery \
	scaffold \
	simple \
	simple-extrap \
	simple-odie \
	simple-rk4tau \
	simple-sbml \
	simple-stoch \
	small-mol \
	tolerance \
	uniMolPtase

DEMO_CONFIG_FILES := $(addprefix $(DEMO)/,$(addsuffix /demo-config.mk, $(DEMO_NAMES)))

include $(DEMO_CONFIG_FILES)
