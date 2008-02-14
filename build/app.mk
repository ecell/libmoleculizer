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

# This differs from module.mk in not having any headers; this might
# have to change.
APP_SOURCE_DIR := $(APP)/$(APP_NAME)
APP_SOURCE_CC_DIR := $(APP_SOURCE_DIR)/cc
APP_SOURCE_CC_FILES := $(addprefix $(APP_SOURCE_CC_DIR)/,$(CC_FILES))
APP_SOURCE_HH_DIR := $(APP_SOURCE_DIR)/hh
APP_SOURCE_HH_FILES := $(addprefix $(APP_SOURCE_HH_DIR)/,$(HH_FILES))

PREEN_LIST += $(APP_SOURCE_DIR)/*~ \
	$(APP_SOURCE_CC_DIR)/*~ \
	$(APP_SOURCE_HH_DIR)/*~

TAGS_LIST := $(APP_SOURCE_HH_FILES) \
	$(TAGS_LIST) \
	$(APP_SOURCE_CC_FILES)

# Dealing with header files.

APP_TARGET_HH_DIR := $(INCLUDE)/$(APP_NAME)
APP_TARGET_HH_FILES := $(addprefix $(APP_TARGET_HH_DIR)/,$(HH_FILES))

$(APP_TARGET_HH_DIR) : | $(INCLUDE)
	mkdir $@

$(APP_TARGET_HH_FILES) : | $(APP_TARGET_HH_DIR)

# This rule's functioning properly depends on the fact that, after a symbolic
# link is created, its modification time is the same as the modification time
# of the file that it points to.  Otherwise this rule would be triggered every
# time a header was modified.  There don't seem to be "static pattern rules"
# for "order only" prerequisites.  This could be an "order only" prerequisite.
$(APP_TARGET_HH_FILES) : $(APP_TARGET_HH_DIR)/%.hh : $(APP_SOURCE_HH_DIR)/%.hh
	ln -s $(PWD)/$< $@

$(INCLUDE)/.headers : | $(APP_TARGET_HH_FILES)


# Dealing with object files and executables.

APP_BUILD_DIR := $(APP_SCRATCH)/$(APP_NAME)

$(APP_BUILD_DIR) : | $(APP_SCRATCH)
	mkdir $@

# Optimized object files and exectuables.
#
# Adds the target app-name/Opt, in particular.  This means that app names
# and module names must be distinct, for now.
EXE_BUILD_DIR := $(APP_BUILD_DIR)/opt

$(EXE_BUILD_DIR) : | $(APP_BUILD_DIR)
	mkdir $@

EXT := Opt
COMPILE_FLAGS := -O2
include $(BUILD)/app-o.mk

# Debugging object files and exectuables.
#
# Adds the target app-name/Dbg, in particular.  This means that app names
# and module names must be distinct, for now.
EXE_BUILD_DIR := $(APP_BUILD_DIR)/dbg

$(EXE_BUILD_DIR) : | $(APP_BUILD_DIR)
	mkdir $@

EXT := Dbg
COMPILE_FLAGS := -g
include $(BUILD)/app-o.mk

# Profiling object files and exectuables.
#
# Adds the target app-name/Dbg, in particular.  This means that app names
# and module names must be distinct, for now.
EXE_BUILD_DIR := $(APP_BUILD_DIR)/prf

$(EXE_BUILD_DIR) : | $(APP_BUILD_DIR)
	mkdir $@

EXT := Prf
COMPILE_FLAGS := -pg
include $(BUILD)/app-o.mk
