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

# Tentatively, this makefile should work in both the source and binary
# releases.

# This copies the part of user's .xmloperator directory that does not
# depend on where the (possibly binary) installation is.
$(USER_XMLOPERATOR) : $(INS)/standalone-xmloperator/.configurable
	rm -rf $@
	cp -r $(INS)/standalone-xmloperator $@

# The directory where the configuration files go.
CONFIG_TGT_DIR := $(USER_XMLOPERATOR)/moleculizer

# This directory is copied from install/xmloperator when user.xmloperator is
# copied, so it doesn't have to be mkdir'ed.
$(CONFIG_TGT_DIR) : | $(USER_XMLOPERATOR)

# Directory containing configuration file stubs, etc.
CONFIG_SRC_DIR := $(INS)/misc

# Absolute paths that go into user configuration files.  CURDIR is
# automatically defined by make, which should always be run from the top-level
# directory.
MOLECULIZER_DIR := $(CURDIR)
XMLO_DIR := $(MOLECULIZER_DIR)/xmloperator_2_3_0
FLAT_DOC_DIR := $(MOLECULIZER_DIR)/doc/static-doc

# Generate the user's shell script to add Moleculizer stuff to PATH and
# LD_LIBRARY_PATH.
$(CONFIG_TGT_DIR)/user.molrc : CSD := $(CONFIG_SRC_DIR)
$(CONFIG_TGT_DIR)/user.molrc : | $(CONFIG_SRC_DIR) $(CONFIG_TGT_DIR)
	echo "#!/bin/bash" > $@
	cat $(CSD)/stub.copyright >> $@
	cat $(CSD)/stub.molrc.begin >> $@
	echo "	MOLECULIZER_DIR=$(MOLECULIZER_DIR)" >> $@
	cat $(CSD)/stub.molrc.end >> $@
	chmod a+rx $@

# Generate the user's shell script to do environment setup for xmloperator.
$(CONFIG_TGT_DIR)/xmlo_rc : CSD := $(CONFIG_SRC_DIR)
$(CONFIG_TGT_DIR)/xmlo_rc : | $(CONFIG_SRC_DIR) $(CONFIG_TGT_DIR)
	echo "#!/bin/bash" > $@
	cat $(CSD)/stub.copyright >> $@
	echo "XMLOPERATOR_HOME=$(XMLO_DIR)" >> $@
	echo "XALAN_LIB=$(XMLO_DIR)" >> $@
	cat $(CSD)/stub.xmlorc >> $@

# Generate thes user's XML configuration file.
# 
# This is used by xslt transformations in the user's xmloperator directory to
# link the HTML translation of an XML input file with the documentation tree
# produced from the schema-doc file that defines the syntax of the input file.
$(CONFIG_TGT_DIR)/mzr-defaults.xml :
	echo -n '<?xml version="1.0" encoding="UTF-8"?>' > $@
	echo -n '<moleculizer-defaults>' >> $@
	echo -n '<server url="' >> $@
	echo -n "$(MOLECULIZER_SERVER)" >> $@
	echo -n '"/>' >> $@
	echo -n '<flat-doc-dir url="file:' >> $@
	echo -n "$(FLAT_DOC_DIR)/" >> $@
	echo -n '"/>' >> $@
	echo -n '</moleculizer-defaults>' >> $@

# Add the user configuration targets to the overall install target.
usr-install : $(CONFIG_TGT_DIR)/user.molrc \
	$(CONFIG_TGT_DIR)/xmlo_rc \
	$(CONFIG_TGT_DIR)/mzr-defaults.xml

PREEN_LIST += $(CONFIG_SRC_DIR)/*~
