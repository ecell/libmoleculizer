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

DOT := $(DOT)/misc

# user.molrc configures user's or server's environment for using moleculizer.
$(DOT)/user.molrc : DT := $(DOT)
$(DOT)/user.molrc : $(DOT)/stub.copyright $(DOT)/stub.molrc.begin $(DOT)/stub.molrc.end $(DOT)/other-things
	echo "#!/bin/bash" > $@
	cat $(DT)/stub.copyright >> $@
	cat $(DT)/stub.molrc.begin >> $@
	echo "	MOLECULIZER_DIR=$(MOLECULIZER_DIR)" >> $@
	echo "	MOLECULIZER_SERVER=$(MOLECULIZER_SERVER)" >> $@
	cat $(DT)/stub.molrc.end >> $@
	chmod a+rx $@

# xmlo_rc configures user's environment for using xmloperator.
$(DOT)/xmlo_rc : DT := $(DOT)
$(DOT)/xmlo_rc : $(DOT)/stub.copyright $(DOT)/stub.xmlorc $(DOT)/other-things
	echo "#!/bin/bash" > $@
	cat $(DT)/stub.copyright >> $@
	echo "XMLOPERATOR_HOME=$(XMLO_DIR)" >> $@
	echo "XALAN_LIB=$(XMLO_DIR)" >> $@
	cat $(DT)/stub.xmlorc >> $@

# Environmental stuff for access from XML.  Mainly conveys server
# URL or flat documentation directory for XSLT transformers to generate
# "live" documentation in xmloperator.
$(DOT)/mzr-defaults.xml : $(DOT)/other-things
	echo -n '<?xml version="1.0" encoding="UTF-8"?>' > $@
	echo -n '<moleculizer-defaults>' >> $@
	echo -n '<server url="' >> $@
	echo -n "$(MOLECULIZER_SERVER)" >> $@
	echo -n '"/>' >> $@
	echo -n '<flat-doc-dir url="file:' >> $@
	echo -n "$(FLAT_DOC_DIR)/" >> $@
	echo -n '"/>' >> $@
	echo -n '</moleculizer-defaults>' >> $@

# Dispatch script for all cgi connected with Moleculizer.  Sets up the
# environment, then invokes mzr-dispatch-task.
$(DOT)/mzr-dispatch : DT := $(DOT)
$(DOT)/mzr-dispatch : $(DOT)/stub.copyright $(DOT)/stub.dispatch $(DOT)/other-things
	echo "#!/bin/bash" > $@
	cat $(DT)/stub.copyright >> $@
	echo ". $(MISC_DIR)/user.molrc" >> $@
	echo "SERVER_HTDOCS_DIR=$(SERVER_HTDOCS_DIR)" >> $@
	echo "SERVER_CGI_DIR=$(SERVER_CGI_DIR)" >> $@
	cat $(DT)/stub.dispatch >> $@
	chmod a+rx $@

# This forces the above scripts to rebuild every time.
$(DOT)/other-things :

# A more selective clean list for cleaning up a "local" installation,
# which allows "power-users" to use their development area
# as their installation.  This allows cleaning of installation-related
# files without complete cleaning.
INSTALL_CLEAN_LIST := $(INSTALL_CLEAN_LIST) \
	$(DOT)/user.molrc \
	$(DOT)/xmlo_rc \
	$(DOT)/mzr-defaults.xml \
	$(DOT)/mzr-dispatch

CLEAN_LIST := $(CLEAN_LIST) \
	$(DOT)/user.molrc \
	$(DOT)/xmlo_rc \
	$(DOT)/mzr-defaults.xml \
	$(DOT)/mzr-dispatch

PREEN_LIST := $(PREEN_LIST) $(DOT)/*~

DOT := $(call dotdot,$(DOT))
