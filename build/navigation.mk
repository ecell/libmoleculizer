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

# Top-level target directories created at build time.
#
# There has to be an "output directory" for XML in order to generate the user's
# .xmloperator directory in a binary distribution, where no source or scratch
# material exists.
LIB := ./lib
BIN := ./bin
DOC := ./doc

# This is a necessary part of the binary release; it has to contain everything
# necessary to make user.xmloperator, aside from makefiles, so much of it will
# just be copied from src/install.  But the schemas and xslt transformations
# can be copied into the skeleton .xmloperator.
# 
# There's no reason not to use it in the same way in the main build, too.
# 
# Now its not clear after all that I need the XML target directory.
INS := ./install
XML := ./xml

# Demos vs. work in progress vs. examples.  Demo->Regression?
DEMO := ./demo
EXAMPLES := ./examples

# The target of install.
USER_XMLOPERATOR := ./user.xmloperator

BUILD := ./build

# Where .o's, .d's, .so's, etc. are generated; same structure as src.
SCRATCH := $(BUILD)/scratch
APP_SCRATCH := $(SCRATCH)/app
UNIT_SCRATCH := $(SCRATCH)/unit
XML_SCRATCH := $(SCRATCH)/xml

# Where everything comes from.
SRC := ./src
INCLUDE := $(SRC)/include
APP := $(SRC)/app
UNIT := $(SRC)/unit
# DOC_SOURCE instead?
DOCS := $(SRC)/doc
INSTALL := $(SRC)/install
# XML_SOURCE instead?
XMLS := $(SRC)/xml
XSL := $(XMLS)/xsl
SCHEMA := $(XMLS)/schema
