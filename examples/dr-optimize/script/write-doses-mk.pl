#!/usr/bin/perl
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

# This version is "doctored" for the dose-response thing to be the top-level
# activity.  Changed DR_OUTPUT_DIR to OUTPUT_DIR.

($doses_file) = @ARGV;

open(DOSES, $doses_file) or
    die("write-doses-mk.pl: Could not open doses file $doses_file.");

# The directories for the respective doses, which appear in the first column
# of the doses file, are all that we need to write the make input, so we skip
# the header line.
$_ = <DOSES>;
($dummy, $dose_var) = split;

while(<DOSES>)
{
    ($dose_dir, $dose) = split;

    print "\nDOSE_DIR := \$(OUTPUT_DIR)/$dose_dir\n";
    print "DOSE_VAR := $dose_var\n";
    print "DOSE := $dose\n";
    print "include \$(BUILD)/dose.mk\n";
}

close(DOSES);
