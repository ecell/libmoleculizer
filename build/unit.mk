###############################################################################
# Nu - a C++ friendly Scheme byte-code compiler.
# Copyright (C) 2004  Walter Lawrence (Larry) Lok.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
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

UNIT_SOURCE_DIR := $(UNIT)/$(UNIT_NAME)
UNIT_SOURCE_CC_DIR := $(UNIT_SOURCE_DIR)/cc
UNIT_SOURCE_CC_FILES := $(addprefix $(UNIT_SOURCE_CC_DIR)/,$(CC_FILES))
UNIT_SOURCE_HH_DIR := $(UNIT_SOURCE_DIR)/hh
UNIT_SOURCE_HH_FILES := $(addprefix $(UNIT_SOURCE_HH_DIR)/,$(HH_FILES))

PREEN_LIST += $(UNIT_SOURCE_DIR)/*~ \
	$(UNIT_SOURCE_CC_DIR)/*~ \
	$(UNIT_SOURCE_HH_DIR)/*~

# Doing the TAGS_LIST like this makes it better for ebrowse, putting
# the header files up front.
TAGS_LIST := $(UNIT_SOURCE_HH_FILES) \
	$(TAGS_LIST) \
	$(UNIT_SOURCE_CC_FILES)

# Dealing with header files.

UNIT_TARGET_HH_DIR := $(INCLUDE)/$(UNIT_NAME)
UNIT_TARGET_HH_FILES := $(addprefix $(UNIT_TARGET_HH_DIR)/,$(HH_FILES))

$(UNIT_TARGET_HH_DIR) : | $(INCLUDE)
	mkdir $@

$(UNIT_TARGET_HH_FILES) : | $(UNIT_TARGET_HH_DIR)

# This rule's functioning properly depends on the fact that, after a symbolic
# link is created, its modification time is the same as the modification time
# of the file that it points to.  Otherwise this rule would be triggered every
# time a header was modified.  There don't seem to be "static pattern rules"
# for "order only" prerequisites.  This could be an "order only" prerequisite.
$(UNIT_TARGET_HH_FILES) : $(UNIT_TARGET_HH_DIR)/%.hh : $(UNIT_SOURCE_HH_DIR)/%.hh
	ln -s $(PWD)/$< $@

$(INCLUDE)/.headers : | $(UNIT_TARGET_HH_FILES)


# Dealing with object files and libraries.

UNIT_BUILD_DIR := $(UNIT_SCRATCH)/$(UNIT_NAME)

$(UNIT_BUILD_DIR) : | $(UNIT_SCRATCH)
	mkdir $@

# Optimized object files and libraries.
#
# Adds the target module-name/Opt, in particular.  App names and module
# names must be distinct, for now.

LIB_BUILD_DIR := $(UNIT_BUILD_DIR)/opt

$(LIB_BUILD_DIR) : | $(UNIT_BUILD_DIR)
	mkdir $@

EXT := Opt
COMPILE_FLAGS := -O2
include $(BUILD)/unit-o.mk

# Debug object files and libraries.
#
# Adds the target module-name/Dbg, in particular.  App names and module
# names must be distinct, for now.

LIB_BUILD_DIR := $(UNIT_BUILD_DIR)/dbg

$(LIB_BUILD_DIR) : | $(UNIT_BUILD_DIR)
	mkdir $@

EXT := Dbg
COMPILE_FLAGS := -g
include $(BUILD)/unit-o.mk

# Profiling object files and libraries.
#
# Adds the target module-name/Prf, in particular.  App names and module
# names must be distinct, for now.

LIB_BUILD_DIR := $(UNIT_BUILD_DIR)/prf

$(LIB_BUILD_DIR) : | $(UNIT_BUILD_DIR)
	mkdir $@

EXT := Prf
COMPILE_FLAGS := -pg
include $(BUILD)/unit-o.mk
