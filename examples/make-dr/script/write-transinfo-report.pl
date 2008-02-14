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

($transinfo_file_stem, $variations_file, $var_suite_output_dir) = @ARGV;

print "write-transinfo-report.pl: transinfo file stem is $transinfo_file_stem.\n";
print "write-transinfo-report.pl: variations file is $variations_file.\n";
print "write-transinfo-report.pl: variation suite outputdirectory is $var_suite_output_dir.\n";

$report_file = "$var_suite_output_dir/tddrs/$transinfo_file_stem--report";
open(REPORT, ">$report_file") or
    die("write-transinfo-report.pl: could not open report file $report_file.");

# Go out to all the variations and get their transinformations, putting
# them into the report and preparing to calculate descriptive statistics.
open(VARIATIONS, $variations_file) or
    die("write-transinfo-report.pl: could not open variations file $variations_file.");
# Skip the header.
<VARIATIONS>;
$total = 0.0;
$tot_sq = 0.0;
$count = 0;
while(<VARIATIONS>)
{
    ($variation_dir) = split;

    $variation_transinfo_file = "$var_suite_output_dir/$variation_dir/tddrs/$transinfo_file_stem--transinfo";
    open(VARIATION_TRANSINFO, $variation_transinfo_file) or
	die("write-transinfo-report.pl: could not open variation transinfo file $variation_transinfo_file.");

    $variation_transinfo = <VARIATION_TRANSINFO>;
    close(VARIATION_TRANSINFO);

    chomp $variation_transinfo;
    print REPORT "$variation_dir\t$variation_transinfo\n";
    $total += $variation_transinfo;
    $tot_sq += ($variation_transinfo * $variation_transinfo);
    $count++;
}
close(VARIATIONS);

# Calculate and print the descriptive statistics in the report.
$mean_transinfo = $total / $count;
$second_moment_transinfo = $tot_sq / $count;
$variance_transinfo = $count * ($second_moment_transinfo - ($mean_transinfo * $mean_transinfo)) / ($count - 1);
$std_dev_transinfo = sqrt($variance_transinfo);
print REPORT "\nmean:\t$mean_transinfo\tstd.dev.:\t$std_dev_transinfo\n";

# Include the transinformation from the mixed population in the report.
$mixed_transinfo_file = "$var_suite_output_dir/tddrs/$transinfo_file_stem--transinfo";
open(MIXED, $mixed_transinfo_file) or
    die("write-transinfo-report.pl: could not open mixed transinfo file $mixed_transinfo_file.");
$mixed_transinfo = <MIXED>;
close(MIXED);
chomp $mixed_transinfo;
print REPORT "\nmixed:\t$mixed_transinfo\n";

close(REPORT);
