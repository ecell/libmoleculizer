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

($variations_file, $output_dir) = @ARGV;

open(VARIATIONS, $variations_file) or
    die("write-variations-mk.pl: could not open variations "
	. "file $variations_file.");

$variations_header = <VARIATIONS>;

while(<VARIATIONS>)
{
    @variation_tokens = split;
    $variation_dir_stem = shift @variation_tokens;
    $variation_dir = "$output_dir/$variation_dir_stem";

    print "DR_OUTPUT_DIR := $variation_dir\n";
    print "include \$(BUILD)/dose-response.mk\n\n";
    print "variations-target : \$(DR_OUTPUT_DIR)/tddr-target\n\n";

    $variation_index++;
}
