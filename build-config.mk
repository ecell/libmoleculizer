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

# Configuration stuff for use at build time.

# xml parsing libraries, libxml2 and libxml++.
XML2_INCLUDE_DIR := /usr/include/libxml2
XML2_LIB := xml2
XMLPP_INCLUDE_DIR := /usr/local/include/libxml++-1.0
XMLPP_LIB_DIR := /usr/local/lib
XMLPP_LIB := xml++-1.0

# Archive library that includes name demangling.
LIBIBERTY_A := /usr/lib/libiberty.a

# Compile and link flags.

CXXFLAGS := -I include -I $(XML2_INCLUDE_DIR) -I $(XMLPP_INCLUDE_DIR) -fPIC -Wall -D_REENTRANT -DPLATFORM_LINUX

SHOB_LINK_FLAGS := -shared -L$(LIB) -L/usr/local/lib -L$(XMLPP_LIB_DIR)

APP_LINK_FLAGS := -L$(LIB) -L/usr/local/lib -L$(XMLPP_LIB_DIR)

