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

# Files to be cleaned away.
CLEAN_LIST := $(SCRATCH) \
	$(BIN) \
	$(LIB) \
	$(XML) \
	$(DOC) \
	$(INS) \
	$(USER_XMLOPERATOR)

# Miscellaneous files (mostly *~'s) to be preened away.
PREEN_LIST := *~ \
	$(BUILD)/*~ \
	$(INCLUDE)

# Now just putting $(BUILD)/unit and $(BUILD)/app on the clean list.
#
# I probably want to generate the lib and bin directories as part of the
# make, and remove them completely in clean.
clean :
	rm -rf $(CLEAN_LIST)

preen : clean
	rm -rf $(PREEN_LIST)
