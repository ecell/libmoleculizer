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

# This program collects the calculated transinfos from the cell types,
# as well as the merged transinfos, and makes a table out of them.

use File::Copy;
use File::Basename;

# Parse the command line.
($values_file_path, $table_output_path) = @ARGV;

# Parse the file names.  

# The job for each line in the values file is executed in a directory whose
# path relative to the VALUES FILE (not the input file) is given in the first
# column of the line.
($values_file_name, $values_file_dir, $values_file_suffix)
    = fileparse($values_file_path);

# Open the table for output.
open(TABLE_HANDLE, ">$table_output_path")
    or die("make-transinfo-table.pl: could not open output file $table_output_path");

# The first line of the tab_values_file gives the property values (like
# VAR1, VAR2, etc.) that are to be substituted with values on remaining lines
# of the tab_values_file.
open(VALUE_HANDLE, $values_file_path);

# Skip the header line.
<VALUE_HANDLE>;

$total_transinfo = 0.0;
$total_sq_transinfo = 0.0;
$transinfo_count = 0;
while(<VALUE_HANDLE>)
{
    @values_line = split;

    $sim_dir_relative_path = shift @values_line;

    $sim_dir = "$values_file_dir/$sim_dir_relative_path";

    $sim_transinfo_file = "$sim_dir/histos/transinfo";

    open(INFO_HANDLE, $sim_transinfo_file) or
	die("make-transinfo-table.pl: could not open transinfo file $sim_transinfo_file");
    $transinfo = <INFO_HANDLE>;
    close(INFO_HANDLE);

    chomp $transinfo;
    print TABLE_HANDLE "$sim_dir_relative_path\t$transinfo\n";
    $total_transinfo += $transinfo;
    $total_sq_transinfo += ($transinfo * $transinfo);
    $transinfo_count++;
}

close(VALUE_HANDLE);

# Put in the average of the transinformations of the cell types.
$average_transinfo = $total_transinfo / $transinfo_count;
$second_moment_transinfo = $total_sq_transinfo / $transinfo_count;
$variance_transinfo
    = $transinfo_count
    * ($second_moment_transinfo - ($average_transinfo
				   * $average_transinfo))
    / ($transinfo_count - 1);
$std_dev_transinfo = sqrt($variance_transinfo);
print TABLE_HANDLE "\n\tMean:\t$average_transinfo\tStd.Dev.: $std_dev_transinfo\n";

# Now put in the transinformation of the mixture.

# This seems unnecessarily hard-wired.
$summary_dir = "transinfo.out/histos";
$summary_transinfo_file = "$summary_dir/transinfo";

open(INFO_HANDLE, "$summary_transinfo_file") or
    die("make-transinfo-table.pl: could not open mixture transinfo file $summary_transinfo_file");

$transinfo = <INFO_HANDLE>;
close(INFO_HANDLE);

chomp $transinfo;
print TABLE_HANDLE "\tMixture:\t$transinfo\n";

$percent_loss = ($average_transinfo - $transinfo) * 100 / $average_transinfo;

print TABLE_HANDLE "\n\tPercent loss: $percent_loss\n";

close(TABLE_HANDLE);

    
    
