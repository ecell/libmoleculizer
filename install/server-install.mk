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

# Get tools for hierarchical approach to make and
# principal landmarks in build tree.
include navigation.mk

# Get server configuration.
include server-config.mk

# Installation configurations: things that depend on where
# the moleculizer ball'o'wax is located (MOLECULIZER_DIR)
MOLECULIZER_DIR := $(shell pwd)
XMLO_DIR := $(MOLECULIZER_DIR)/xmloperator_2_3_0
FLAT_DOC_DIR := $(MOLECULIZER_DIR)/doc/static-doc
MISC_DIR := $(MOLECULIZER_DIR)/install/misc

# Install moleculizer into web server's htdocs and cgi directories.
# The web server must be configured to follow symbolic links.
# This needs to be done as root, or another user with permission
# to put things in htdocs and cgi-bin.
#
# This automatic use of ldconfig will have to be changed, in all probability.
server-install : server/install-target $(INSTALL)/server-xmloperator/install-target $(SERVER_CGI_DIR)/mzr-dispatch $(SERVER_HTDOCS_DIR)/moleculizer
	rm -rf user.xmloperator
	cp -a $(INSTALL)/server-xmloperator user.xmloperator
	ldconfig $(MOLECULIZER_DIR)/lib

# Remove stuff installed in www's htdocs and cgi directories.
server-remove :
	rm -rf $(SERVER_HTDOCS_DIR)/moleculizer
	rm -f $(SERVER_CGI_DIR)/mzr-dispatch

# Install cgi dispatch script.
$(SERVER_CGI_DIR)/mzr-dispatch : $(INSTALL)/misc/mzr-dispatch other-things
	cp $< $@

# Copy the server directory to www's htdocs directory.
$(SERVER_HTDOCS_DIR)/moleculizer : server other-things
	rm -rf $@
	cp -r server $@
	chmod -R a+rwx $@

INSTALL_CLEAN_LIST := user.xmloperator
install-clean :
	rm -rf $(INSTALL_CLEAN_LIST)

other-things :

include $(INSTALL)/td.mk \
	$(SERVER)/td.mk \
	demo/td.mk
