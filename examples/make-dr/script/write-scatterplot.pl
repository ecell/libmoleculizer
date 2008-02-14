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

# Emits DORA scores and transinformations from all the variations into a
# single file for later plotting as a scatterplot of DORA score against
# transinformation.  

# There is a different transinformation for each response, but only one DORA
# score, which may involve all the responses.  Thus, I end up make a
# scatterplot for each "response."

($transinfo_file_stem, $variations_file, $var_suite_output_dir) = @ARGV;

print "write-scatterplot.pl: transinfo file stem is $transinfo_file_stem.\n";
print "write-scatterplot.pl: variations file is $variations_file.\n";
print "write-scatterplot.pl: variation suite output directory is $var_suite_output_dir.\n";

$scatterplot_file = "$var_suite_output_dir/tddrs/$transinfo_file_stem--scatterplot";
open(SCATTERPLOT, ">$scatterplot_file") or
    die("write-scatterplot.pl: could not open scatterplot file $scatterplot_file.");

# Go out to all the variations and get their transinformations, putting
# them into the report and preparing to calculate descriptive statistics.
open(VARIATIONS, $variations_file) or
    die("write-scatterplot.pl: could not open variations file $variations_file.");
# Skip the header.
<VARIATIONS>;
while(<VARIATIONS>)
{
    ($variation_dir) = split;

    $variation_transinfo_file = "$var_suite_output_dir/$variation_dir/tddrs/$transinfo_file_stem--transinfo";
    open(VARIATION_TRANSINFO, $variation_transinfo_file) or
	die("write-scatterplot.pl: could not open variation transinfo file $variation_transinfo_file.");

    $variation_transinfo = <VARIATION_TRANSINFO>;
    close(VARIATION_TRANSINFO);
    chomp $variation_transinfo;

    $variation_score_file = "$var_suite_output_dir/$variation_dir/tddrs/dora-score";
    open(VARIATION_SCORE, $variation_score_file) or
	die("write-scatterplot.pl: could not open variation score file $variation_score_file.");
    $variation_score = <VARIATION_SCORE>;
    close(VARIATION_SCORE);
    chomp $variation_score;

    print SCATTERPLOT "$variation_score\t$variation_transinfo\n";
}
close(VARIATIONS);
close(SCATTERPLOT);
