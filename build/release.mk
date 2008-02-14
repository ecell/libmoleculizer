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

# Rules to generate binary and source releases.

BINARY_RELEASE_NAME := $(PROJECT)-$(MAJOR_VERSION).$(MINOR_VERSION)-binary

BINARY_RELEASE_TARBALL := $(BINARY_RELEASE_NAME).tgz

TMP_BINARY_RELEASE := /tmp/$(BINARY_RELEASE_NAME)

# After everything is done to prepare the binary release off in /tmp, it is
# made into a tar file here.
$(BINARY_RELEASE_TARBALL) : $(TMP_BINARY_RELEASE)/.prepared
	cd /tmp \
	&& rm -f $@ \
	&& tar -czf $@ $(BINARY_RELEASE_NAME)
	cp /tmp/$@ .

# This target also depends on getting the user-xmloperator stuff together.
.binary-release-prepared : opt xml usr-doc

# Note that .prepared is a phony target, forcing this to rebuild every
# time.
# 
# Note that this makes a rather funky use of the fact that the directory
# trees are the same.
# 
# Note that this project should be 'svn export'-ed before doing the release.
# It could be built into this script, but then the source release wouldn't be
# able to generate a binary release.  Besides, it appears that svn export does
# funny things with modification times, apparently in order to accomodate
# programmers in different time zones.
$(TMP_BINARY_RELEASE)/.prepared :
	rm -rf $(@D)
	cp -a . $(@D)
	cd $(@D) \
	&& $(MAKE) -f newMakefile .binary-release-prepared \
	&& rm -rf src \
	&& rm -rf $(SCRATCH) \
	&& cp build/binary-release-makefile ./newMakefile

binary-release : $(BINARY_RELEASE_TARBALL)


######################

SOURCE_RELEASE_NAME := $(PROJECT)-$(MAJOR_VERSION).$(MINOR_VERSION)-source

SOURCE_RELEASE_TARBALL := $(SOURCE_RELEASE_NAME).tgz

TMP_SOURCE_RELEASE := /tmp/$(SOURCE_RELEASE_NAME)

$(SOURCE_RELEASE_TARBALL) : $(TMP_SOURCE_RELEASE)/.prepared
	cd /tmp \
	&& rm -f $@ \
	&& tar -czf $@ $(SOURCE_RELEASE_NAME)
	cp /tmp/$@ .

# Note that this makes a rather funky use of the fact that the directory
# trees are the same.
# 
# Note that this project should be preened and 'svn export'-ed before doing
# the release.  Note that only the scratch directory is removed before
# tarring.
 $(TMP_SOURCE_RELEASE)/.prepared :
	rm -rf $(@D)
	cp -a . $(@D)
	cd $(@D) \
	&& rm -rf $(SCRATCH)

source-release : $(SOURCE_RELEASE_TARBALL)