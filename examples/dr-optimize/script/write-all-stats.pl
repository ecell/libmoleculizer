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

($histo_dir, $doses_file, $all_stats_file) = @ARGV;

open(OUTPUT, ">$all_stats_file") or
    die("write-all-stats.pl: could not open output file $all_stats_file.");

open(DOSES, $doses_file) or
    die("write-all-stats.pl: could not open doses file $doses_file.");

# Skip the header line of the doses file.
<DOSES>;

while(<DOSES>)
{
    ($dose_dir, $dose) = split;

    $data_file = "$histo_dir/$dose.data";
    open(DATA, $data_file) or
	die("write-all-stats.pl: could not open data file $data_file.");

    $total = 0.0;
    $totSq = 0.0;
    $count = 0;
    while($value = <DATA>)
    {
	chomp $value;

	$total += $value;
	$totSq += ($value * $value);
	$count++;
    }
    close(DATA);

    $mean = $total / $count;
    $second_moment = $totSq / $count;
    $sampleVariance
	= $count * ($second_moment - ($mean * $mean)) / ($count - 1);
    $standard_deviation = sqrt($sampleVariance);
    if(0.0 < $mean)
    {
	$relative_error = $standard_deviation / $mean;
    }
    else
    {
	# This is meaningless in any case, so if we have any use for it
	# one will have to come back to this point.
	#
	# Relative error was included at Andrew's request.
	$relative_error = 0.0;
    }

    # Just emitting mean and standard deviation might make more sense here.
    print OUTPUT "$dose\t$mean\t$sampleVariance\t$relative_error\n";
}

close(DOSES);

close(OUTPUT);
