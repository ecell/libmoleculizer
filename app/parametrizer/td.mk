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

APP_NAME := parametrizer

# Units compiled by this project that are required as shared objects
# (or archives for profiling.)
REQUIRED_UNITS := domUtils \
	sampleDist \
	mzr \
	mol \
	plex \
	dimer \
	stoch \
	bndKinase

# REQUIRED_UNITS := domUtils \
# 	sampleDist \
# 	mzr \
# 	mol \
# 	plex \
# 	dimer \
# 	stoch \
# 	gpa \
# 	nucEx \
# 	modKinase \
# 	scaffold \
# 	bndKinase

# Libraries not compiled by this project, required as shared objects.
EXTRA_LIBS :=

# List of archive libraries for static linking.  This list may (usually does)
# contain repetitions, and order is important.  Normally, every element of
# this list should appear in 'REQUIRED_UNITS' above, which is used to
# construct dependencies on the archives. This variable is used to make the
# profiling version of the program.
STATIC_ARCHIVE_LINK_LIST := mzr \
	bndKinase \
	scaffold \
	nucEx \
	modKinase \
	gpa \
	plex \
	dimer \
	mol \
	stoch \
	mzr \
	domUtils \
	sampleDist

# List of libraries that have to be linked dynamically into statically-built
# executable.  This is normally the union of all the "EXTRA_LIBS" for all
# the "REQUIRED_UNITS" above.  In the dynamic-linking regime, these libraries
# are linked into the unit .so's.  In the static-linking regime, the only
# place for them to go in is when the executable is being linked.
STATIC_EXTRAS := xml++-1.0 \
	xml2

DOT := $(DOT)/$(APP_NAME)

include $(DOT)/cc/td.mk
include $(DOT)/dbg-o/td.mk
include $(DOT)/opt-o/td.mk
include $(DOT)/prf-o/td.mk

PREEN_LIST := $(PREEN_LIST) $(DOT)/*~

DOT := $(call dotdot,$(DOT))
