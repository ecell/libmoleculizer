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

# Generation of user's .xmloperator directory, source release, and binary
# release.

# This target copies the entire moleculizer ball'o'wax to /tmp,
# then removes soureces, etc. from the copy, then moves the copy
# into this directory as the standalone binary delivery.
#
# Note that this doesn't do anything with svn.  One should use
# svn export to get rid of .svn files, then make this target.
STANDALONE_BINARY_DIR := $(PROJECT)-standalone-binary-$(MAJOR_VERSION).$(MINOR_VERSION)
STANDALONE_STAGING_DIR := /tmp/$(STANDALONE_BINARY_DIR)
standalone-binary : opt
	rm -rf $(STANDALONE_BINARY_DIR)
	rm -rf $(STANDALONE_STAGING_DIR)
	cp -a . $(STANDALONE_STAGING_DIR)
	$(MAKE) -C $(STANDALONE_STAGING_DIR) standalone-clean
	cp -a $(STANDALONE_STAGING_DIR) $(STANDALONE_BINARY_DIR)
	rm -rf $(STANDALONE_STAGING_DIR)

# Removes what is not needed in a standalone binary release. (Some XSL files
# are needed to do transformations in the demos.) Then swaps in new makefile.
standalone-clean :
	rm -rf app include unit
	rm -rf server session-makefiles
	mv $(INSTALL)/standalone-install.mk makefile

# Generates preened source delivery.
#
# Note that this doesn't do anything with svn.  One should use
# svn export to get rid of .svn files, then make this target.
SOURCE_RELEASE_DIR := $(PROJECT)-source-$(MAJOR_VERSION).$(MINOR_VERSION)
SOURCE_STAGING_DIR := /tmp/$(SOURCE_RELEASE_DIR)
source-release :
	rm -rf $(SOURCE_RELEASE_DIR)
	rm -rf $(SOURCE_STAGING_DIR)
	cp -a . $(SOURCE_STAGING_DIR)
	$(MAKE) -C $(SOURCE_STAGING_DIR) preen
	cp -a $(SOURCE_STAGING_DIR) $(SOURCE_RELEASE_DIR)
	rm -rf $(SOURCE_STAGING_DIR)

# For using the development area as an installation.
local-install : install/target
	$(MAKE) -f $(INSTALL)/standalone-install.mk standalone-install
