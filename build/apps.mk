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

APP_NAMES := continuator \
	cptzr \
	moleculizer \
	odie \
	parametrizer \
	rk4tau \
	xmlpretty

PREEN_LIST += $(APP)/*~ \
	$(APP)/scripts/*~

# Target directory for executables.
$(BIN) :
	mkdir $@

# Scratch directory for builds of all apps.
$(APP_SCRATCH) : | $(SCRATCH)
	mkdir $@

include $(addprefix $(APP)/,$(addsuffix /app-config.mk,$(APP_NAMES)))

SCRIPTS := altSeedString \
	concat.pl \
	getSeedString \
	gpFunctions.pm \
	mzr-dispatch-task \
	plot-cnc-files \
	plot-dmp-files \
	plot-dmp-files-as-cnc \
	plt \
	ready-targets.bash \
	state2odie \
	state2rk4tau \
	state2sbml \
	tddr-int-setup.pl \
	tddr-int-setup-odie.pl \
	tddr-extract-timepoint.pl \
	tddr-float-setup.pl \
	tddr-float-setup-odie.pl \
	tddr-plot-timepoint.sh \
	tddr-report.pl \
	vary-param

SOURCE_SCRIPT_DIR := $(APP)/scripts
SOURCE_SCRIPTS := $(addprefix $(SOURCE_SCRIPT_DIR)/,$(SCRIPTS))
TARGET_SCRIPTS := $(addprefix $(BIN)/,$(SCRIPTS))

$(TARGET_SCRIPTS) : $(BIN)/% : $(SOURCE_SCRIPT_DIR)/% | $(BIN)
	cp $< $@
	chmod a+x $@

# Add the applications and the scripts to the basic targets.
opt : $(addsuffix /Opt/dynamic,$(APP_NAMES)) $(TARGET_SCRIPTS)

dbg : $(addsuffix /Dbg/dynamic,$(APP_NAMES)) $(TARGET_SCRIPTS)
