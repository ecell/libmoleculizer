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

# The default target and some documentation of other primary
# targets.

default : opt xml-target usr-doc

# Produces dynamically-linked, optimized executables.
opt : 

# Produces dynamically-linked, debugging executables.
dbg :

# Produces statically-linked executables for profiling,
# performance analysis.
prf :

# Produces various kinds of schema files from schema-doc files; copies xslt
# transformations.
xml-target :

# Produces trees of HTML documentation of XML file formats from
# schema-doc files.
usr-doc :

# Produces documentation of source code, using Doxygen.
src-doc :

# Generates the user.xmloperator directory, which each user installs in her
# own home directory in order to use Moleculizer, xmloperator, etc.
usr-install :

# Generates a binary release.
# 
# First, the code is exported from Subversion to remove all the .svn
# directories from the tree, followed by the default build.  Then src and
# build/scratch are deleted.  Then the top-level makefile is replaced with one
# that does user-install, demos, and examples, but that doesn't look for
# src.
# 
# Then, the whole thing is made into a tarball with the right name.
binary-release :

# Generates a source release.
#
# First, the code is exported from Subversion to remove all the .svn
# directories from the tree.
# 
# Then, the whole thing is made into a tarball with the right name.
source-release : 

# This triggers generation and loading of all the dependency files.
test-target :
	@echo Dummy target compiled!

######################################################

# Include landmarks in the build tree.
include build/navigation.mk

# Include build configuration, such as library locations and
# compiler flags, configuration for different platforms, etc.
include build/build-config.mk

# Include rules for removing products of build; clean and preen.
include build/clean.mk

# Include rules for generating tags files and browse files.  This file list
# is also used for generating html source-code documentation using doxygen.
include build/tags.mk

# Include rules for assembling src/include during the build and removing
# it at clean time.
include build/headers.mk

# Include rules for building, cleaning, etc. shared libraries (units).
include build/units.mk

# Include rules for building, cleaning, etc. apps.
include build/apps.mk

# Include rules for processing schemas and other xml files.
include build/xml.mk

# Include rules for processing documentation.
include build/doc.mk

# Add rules to generate user's .xmloperator directory, which performs
# all user's environment setup for moleculizer, xmloperator, etc.
include build/install.mk

# Add rules to generate binary release.
include build/release.mk

# Add rules to run the demos.
include build/demos.mk

# Add rules to run the examples.
include build/examples.mk

# This arranges that all the headers in $(INCLUDE) are up-to-date with respect
# to their modules before dependency checking begins.  The headers must be in
# place so that g++ can find them in its dependency-generation pass.
$(MAKEFILES) : $(INCLUDE)/.headers

# Include dependency files.
-include $(MAKEFILES)
