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

# Names of local executable targets.
LOCAL_EXE_STEM := $(DOT)/$(APP_NAME)$(EXT)
DYNAMIC_LOCAL_EXE := $(LOCAL_EXE_STEM)Dynamic
STATIC_LOCAL_EXE := $(LOCAL_EXE_STEM)Static

# The executable in bin that is overwritten by the "static" and "dynamic"
# targets.
REMOTE_EXE := $(BIN)/$(APP_NAME)

# Copies local static executable to bin directory with generic name.
$(DOT)/static : RMT_EXE := $(REMOTE_EXE)
$(DOT)/static : $(STATIC_LOCAL_EXE)
	cp $< $(RMT_EXE)

# Copies local static executable to bin directory with generic name.
$(DOT)/dynamic : RMT_EXE := $(REMOTE_EXE)
$(DOT)/dynamic : $(DYNAMIC_LOCAL_EXE)
	cp $< $(RMT_EXE)

OBJECTS := $(addprefix $(DOT)/,$(SOURCES:.cc=.o))


########## Local statically-linked executable.

# Input archive libraries built in this project, for use as a prerequitie.
REQUIRED_ARCHIVES := $(addprefix $(LIB)/lib,$(addsuffix $(EXT).a,$(REQUIRED_UNITS)))

# Input archive library link list, which may have necessary repititions, for
# use in statically linking this app.  Some of these are built in this
# project, some are externally provided.
LNK_MODS := $(addprefix $(LIB)/lib,$(addsuffix $(EXT).a,$(STATIC_ARCHIVE_LINK_LIST)))

# Dynamic libraries that must be linked into the static executable.
# Typically, this is external stuff not supplied in archive form.
LNK_SYS := $(addprefix -l,$(EXTRA_LIBS))

# Statically links local executable.
$(STATIC_LOCAL_EXE) : LNKS := $(OBJECTS) $(LNK_MODS) $(EXTRA_ARCHIVES) $(LNK_SYS)
$(STATIC_LOCAL_EXE) : $(REQUIRED_ARCHIVES) $(OBJECTS) 
	g++ $(APP_LINK_FLAGS) -o $@ $(LNKS)


########## Local dynamically-linked executable.

# Input dynamic libraries built in this project, for use as a prerequisite.
REQUIRED_SOS := $(addprefix $(LIB)/lib,$(addsuffix $(EXT).so,$(REQUIRED_UNITS)))
# Input dynamic libraries with -l prefix, for use in linking.
LNK_MODS := $(addprefix -l,$(addsuffix $(EXT),$(REQUIRED_UNITS)))

# External dynamic libraries.
LNK_SYS  := $(addprefix -l,$(EXTRA_LIBS))

$(DYNAMIC_LOCAL_EXE) : LNKS :=  $(LNK_SYS) $(LNK_MODS) $(OBJECTS) $(EXTRA_ARCHIVES)
$(DYNAMIC_LOCAL_EXE) : $(OBJECTS) $(REQUIRED_SOS)
	g++ $(APP_LINK_FLAGS) -o $@ $(LNKS)


########## Object files and their dependiencies.

# Local C++ compiler flags, all put together.
COMPILE_LOCAL_FLAGS := $(COMPILE_FLAGS) $(COMPILE_DEFS)

# Compiles local objects
$(OBJECTS) : CLF := $(COMPILE_LOCAL_FLAGS)
$(OBJECTS) : $(DOT)/%.o : $(CC_DIR)/%.cc
	$(CXX) -c $(CXXFLAGS) $(CLF) -o $@ $<

# Generates dependencies for this application's objects.
APP_MAKEFILES := $(addprefix $(DOT)/,$(SOURCES:.cc=.d))
$(APP_MAKEFILES) : CLF := $(COMPILE_LOCAL_FLAGS)
$(APP_MAKEFILES) : $(DOT)/%.d : $(CC_DIR)/%.cc
	g++ -MM -MT '$@ $(@:.d=.o)' $(CXXFLAGS)  $(CLF) $< > $@

MAKEFILES := $(MAKEFILES) $(APP_MAKEFILES)


########## Housekeeping.

CLEAN_LIST := $(CLEAN_LIST) \
	$(OBJECTS) \
	$(STATIC_LOCAL_EXE) \
	$(DYNAMIC_LOCAL_EXE) \
	$(REMOTE_EXE)

PREEN_LIST := $(PREEN_LIST) $(DOT)/*~

