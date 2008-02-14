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

# Name of the application as it will be invoked; e.g. moleculizer.
APP_NAME := xmlpretty

# The source files to be compiled, including the .cc extension.  These should
# all be in the cc directory.
CC_FILES := xmlpretty.cc

# Header files for this app.  These go into a subdirectory of include
# with the same name as the app.
HH_FILES :=

# Units compiled by this project that are required for this executable.
# The archive form is used for static linking, and the shared object form
# is used for dynamic linking.
REQUIRED_UNITS := utl

# Shared objects not compiled by this project, linked into both static and
# dynamic versions of the app.
EXTRA_LIBS := $(XML2_LIB) \
	$(XMLPP_LIB)

# List of unit archives for static linking.  This list may (usually does)
# contain repetitions, and order is important.  Normally, every element of
# this list should appear in 'REQUIRED_UNITS' above, which is used to
# construct dependencies on the archives.
STATIC_ARCHIVE_LINK_LIST := domUtils

# Archive libraries, not compiled by this project, linked into both static
# and dynamic versions of this app.
EXTRA_ARCHIVES :=

# Special compiler flags and definitions to be used in compiling c++ files for
# this app. (Rarely set.)
COMPILE_DEFS :=

include $(BUILD)/app.mk
