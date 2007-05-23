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

PLATFORM := PLATFORM_SUSE_LINUX

ifeq ($(PLATFORM), PLATFORM_SUSE_LINUX)

# Demangler headers from libcwd.
DEMANGLE_INCLUDE_DIR := /usr/include/libcwd/

# Gnome's C-based XML parsing library.
XML2_INCLUDE_DIR := /usr/include/libxml2
XML2_LIB_DIR :=
XML2_LIB := xml2

# C++ wrapper library for XML2.
XMLPP_INCLUDE_DIR := /usr/include/libxml++-1.0
XMLPP_LIB_DIR := /usr/lib
XMLPP_LIB := xml++-1.0


GLIBMM_COMPILE_FLAGS := $(shell pkg-config --cflags glibmm-2.4)
# Note that these come with the -L attached.  An alternative would
# be to remove the -L's using make's subst or patsubst function.


GLIBMM_LIB_DIRS := $(shell pkg-config --libs-only-L glibmm-2.4)
# Note that these come with the -l attached.  An alternative would
# be to remove the -l's using make's subst or patsubst function.


GLIBMM_LIBS := $(shell pkg-config --libs-only-l glibmm-2.4)

# Gnu scientific library.
GSL_INCLUDE_DIR :=
GSL_LIB_DIR :=
GSL_LIB := gsl

# GSL linear algebra.
GSL_CBLAS_INCLUDE_DIR :=
GSL_CBLAS_LIB_DIR :=
GSL_CBLAS_LIB := gslcblas

# Archive library that includes name demangling.
# LIBIBERTY_A := /usr/lib/libiberty.a
LIBIBERTY_A := 

endif # PLATFORM_SUSE_LINUX

#################

PROJECT := moleculizer
MAJOR_VERSION := 1
MINOR_VERSION := 0.9

#################

# Compile flags.
INCLUDE_DIRS := $(INCLUDE) \
	$(DEMANGLE_INCLUDE_DIR) \
	$(XML2_INCLUDE_DIR) \
	$(XMLPP_INCLUDE_DIR) \
	$(GSL_INCLUDE_DIR) \
	$(GSL_CBLAS_INCLUDE_DIR)

INCLUDE_DIR_FLAGS := $(addprefix -I ,$(INCLUDE_DIRS))

CXXFLAGS := $(INCLUDE_DIR_FLAGS) -fPIC -Wall -D_REENTRANT -D$(PLATFORM)

#################

# Link flags.
LINK_DIRS := $(LIB) \
	$(XML2_LIB_DIR) \
	$(XMLPP_LIB_DIR) \
	$(GSL_LIB_DIR) \
	$(GSL_CBLAS_LIB_DIR)

LINK_DIR_FLAGS := $(addprefix -L ,$(LINK_DIRS))

SHOB_LINK_FLAGS := $(LINK_DIR_FLAGS) -shared 

APP_LINK_FLAGS := $(LINK_DIR_FLAGS)

# Set java class path
export CLASSPATH=$CLASSPATH:xmloperator_2_3_0/xalan.jar