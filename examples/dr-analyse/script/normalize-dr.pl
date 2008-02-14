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

use POSIX;

# Rescale the dose-response file so that responses are all on a [0,1] scale.

($all_stats_file, $dose_response_file) = @ARGV;

open(STATS, $all_stats_file) or
    die("normalize.pl: could not open all stats file $all_stats_file.");

$max_response = POSIX::DBL_MIN;
$min_response = POSIX::DBL_MAX;
while(<STATS>)
{
    ($dose, $mean_response, $response_variance) = split;

    push @doses, $dose;
    push @responses, $mean_response;

    if($max_response < $mean_response)
    {
	$max_response = $mean_response;
    }
    if($min_response > $mean_response)
    {
	$min_response = $mean_response;
    }
}
close(STATS);

$scale_factor = $max_response - $min_response;

open(DR, ">$dose_response_file") or
    die("normalize.pl: could not open dose-response file $dose_response_file.");
while(@doses)
{
    $dose = shift @doses;
    $response = shift @responses;
    $normalized_response = ($response - $min_response) / $scale_factor;

    print DR "$dose\t$normalized_response\n";
}

close(DR);


    

