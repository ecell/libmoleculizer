###############################################################################
# BNGMZRConvertor - A utility program for converting bngl input files to mzr
#		    input files.
# Copyright (C) 2007 Nathan Addy
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
#   Nathan Addy, Scientific Programmer	Voice: 510-981-8748
#   The Molecular Sciences Institute    Email: addy@molsci.org  
#   2168 Shattuck Ave.                  
#   Berkeley, CA 94704
###############################################################################



COMPILED_PYTHON := $(shell find . -name "*pyc" | grep -v "\.svn")
EMACS_JUNK := $(shell find . -name "*~" -o -name "\#*\#" | grep -v "\.svn")
CLEAN_FILES := $(COMPILED_PYTHON) \
	$(EMACS_JUNK)


default:
	

clean:
	rm $(CLEAN_FILES)