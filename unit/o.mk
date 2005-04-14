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

OBJECTS := $(addprefix $(DOT)/,$(SOURCES:.cc=.o))

LINKER_NAME := lib$(UNIT_NAME)$(EXT).so
SONAME := $(LINKER_NAME).$(MAJOR_VERSION)
LIB_NAME := $(SONAME).$(MINOR_VERSION)
ARCHIVE_NAME := lib$(UNIT_NAME)$(EXT).a

# Paths to required unit libraries (libraries made by this build
# that this unit must link with) for checking that they are up-to-date.
REQUIRED_UNITS_LIBS := $(addprefix $(LIB)/lib,$(addsuffix $(EXT).so,$(REQUIRED_UNITS)))

# Name of units for required unit libraries (libraries made by this build
# that this unit must link with) for use by ld.
REQUIRED_UNITS_LD := $(addprefix -l,$(addsuffix $(EXT),$(REQUIRED_UNITS)))

# Names of "external" libraries (not made by this build) for use by ld.
EXTRA_LIBS_LD := $(addprefix -l,$(EXTRA_LIBS))

# The linker's soname flag for this unit
SONAME_LD := -Wl,-soname,$(SONAME)

# Most link flags, all put together.
LNK_LD := $(SONAME_LD) $(EXTRA_LIBS_LD) $(REQUIRED_UNITS_LD) $(OBJECTS) $(EXTRA_ARCHIVES)

# This is intended for dynamic linking units, which need not be
# prerequisites of any executable.
$(DOT)/dynamic : $(LIB)/$(LINKER_NAME)

# Generates link to unit shared object for use by ld.
$(LIB)/$(LINKER_NAME) : $(LIB)/$(SONAME)
	ln -sf $(<F) $@

# Generates link to unit shared object for use by ld.so.
$(LIB)/$(SONAME) : $(LIB)/$(LIB_NAME)
	ln -sf $(<F) $@

# Links unit shared object.
$(LIB)/$(LIB_NAME) : LNKS := $(LNK_LD)
$(LIB)/$(LIB_NAME) : $(OBJECTS) $(REQUIRED_UNITS_LIBS)
	g++ $(SHOB_LINK_FLAGS) -o $@ $(LNKS)

# Generates archive of unit objects.
$(LIB)/$(ARCHIVE_NAME) : $(OBJECTS)
	ar rc $@ $^
	ranlib $@

# Local C++ compiler flags, all put together.
COMPILE_LOCAL_FLAGS := $(COMPILE_FLAGS) $(COMPILE_DEFS)

# Compiles the objects for this unit.
$(OBJECTS) : CLF := $(COMPILE_LOCAL_FLAGS)
$(OBJECTS) : $(DOT)/%.o : $(DOT_DOT)/cc/%.cc
	$(CXX) -c $(CXXFLAGS) $(CLF) $(DEFS) -o $@ $<

# Generates dependencies for objects in this unit.
UNIT_MAKEFILES := $(addprefix $(DOT)/,$(SOURCES:.cc=.d))
$(UNIT_MAKEFILES) : CLF := $(COMPILE_LOCAL_FLAGS)
$(UNIT_MAKEFILES) : $(DOT)/%.d : $(DOT_DOT)/cc/%.cc
	g++ -MM -MT '$@ $(@:.d=.o)' $(CXXFLAGS)  $(CLF) $< > $@

# 	g++ -MM -MT $(@:.d=.o) $(CXXFLAGS)  $(CLF) $< > $@.$$$$; \
# 	sed 's/^\(.*\)\.o\s*:\s*/\1.d \1.o : /' < $@.$$$$ > $@; \
# 	rm -f $@.$$$$

MAKEFILES := $(MAKEFILES) $(UNIT_MAKEFILES)

CLEAN_LIST := $(CLEAN_LIST) \
	$(OBJECTS) \
	$(LIB)/$(SONAME) \
	$(LIB)/$(LINKER_NAME) \
	$(LIB)/$(LIB_NAME) \
	$(LIB)/$(ARCHIVE_NAME)

PREEN_LIST := $(PREEN_LIST) $(DOT)/*~

