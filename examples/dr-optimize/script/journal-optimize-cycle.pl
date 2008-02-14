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

# This script writes an entry in 'moleculizer-variations' for each optimization
# step.  The file moleculizer-variations can then be used to control a standard
# suite of dose/response sweeps, similar to the bottom two tiers of my usual
# three-tiered sweep of reaction rates, populations, and doses.

# The index of the current optimization step.
$ndx = shift(@ARGV);
# The "size," an error estimate, of the current optimization step.
$current_size = shift(@ARGV);
# The value being optimized, usually transinformation.
$current_value = shift(@ARGV);
# All the parameters.
@parameters = @ARGV;

# This file name is just hardwired in. This was an unexpected advantage of
# using a separate script for journaling.
$variations_file_name = "moleculizer-variations";

# Name of the directory where the "recapitulating" simulation for this
# optimization step will be done.  This recapitulating simulation is supposed
# to do full analysis, dose/response curves, DORA scores, etc.
$variation_dir = "opt-$ndx";

open(VARS, ">>$variations_file_name") or
    die("Could not open variations file $variations_file_name for appending.");
print VARS "$variation_dir\t";
while(@parameters)
{
    # This variations file is for use by the "analyser" which uses plain
    # rates, instead of logarithmic rates.
    $log_rate = shift(@parameters);
    $rate = exp($log_rate);
    # I get too many digits here sometimes....printf?
    print VARS "$rate";

    if(@parameters)
    {
	print VARS "\t";
    }
    else
    {
	print VARS "\n";
    }
}
close(VARS);

