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

DOT_DOT := $(DOT)
DOT := $(DOT)/opt-o

EXT := Opt
COMPILE_FLAGS := -O2

# Copies local executable to bin with generic name.
LOCAL_EXE := $(DOT)/$(APP_NAME)$(EXT)

$(DOT)/target : BIN_EXE := $(BIN)/$(APP_NAME)
$(DOT)/target : $(LOCAL_EXE)
	cp $< $(BIN_EXE)

# Links local executable
OBJECTS := $(addprefix $(DOT)/,$(SOURCES:.cc=.o))
REQUIRED_SOS := $(addprefix $(LIB)/lib,$(addsuffix $(EXT).so,$(REQUIRED_UNITS)))
LNK_MODS := $(addprefix -l,$(addsuffix $(EXT),$(REQUIRED_UNITS)))
LNK_SYS  := $(addprefix -l,$(EXTRA_LIBS))

$(LOCAL_EXE) : LNKS :=  $(LNK_SYS) $(LNK_MODS) $(OBJECTS)
$(LOCAL_EXE) : $(OBJECTS) $(REQUIRED_SOS)
	g++ $(APP_LINK_FLAGS) -o $@ $(LNKS)

# Definitions for inclusion in C++ files.
COMPILE_DEFS :=

# Local C++ compiler flags, all put together.
COMPILE_LOCAL_FLAGS := $(COMPILE_FLAGS) $(COMPILE_DEFS)

# Compiles local objects
$(OBJECTS) : CLF := $(COMPILE_LOCAL_FLAGS)
$(OBJECTS) : $(DOT)/%.o : $(DOT_DOT)/cc/%.cc
	$(CXX) -c $(CXXFLAGS) $(CLF) -o $@ $<

# Generates dependencies for this application's objects.
APP_MAKEFILES := $(addprefix $(DOT)/,$(SOURCES:.cc=.d))
$(APP_MAKEFILES) : CLF := $(COMPILE_LOCAL_FLAGS)
$(APP_MAKEFILES) : $(DOT)/%.d : $(DOT_DOT)/cc/%.cc
	g++ -MM -MT $(@:.d=.o) $(CXXFLAGS) $(CLF) $< > $@.$$$$; \
	sed 's/^\(.*\)\.o\s*:\s*/\1.d \1.o : /' < $@.$$$$ > $@; \
	rm -f $@.$$$$

MAKEFILES := $(MAKEFILES) $(APP_MAKEFILES)

CLEAN_LIST := $(CLEAN_LIST) \
	$(addprefix $(DOT)/,$(SOURCES:.cc=.o)) \
	$(LOCAL_EXE) \
	$(BIN)/$(APP_NAME)

PREEN_LIST := $(PREEN_LIST) $(DOT)/*~

DOT := $(call dotdot,$(DOT))
