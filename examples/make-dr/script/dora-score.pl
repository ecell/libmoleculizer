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

# The idea of this script, beyond what I wrote earlier, is to allow the
# possibility that all the responses (really, their Hill function parameters)
# can be taken into account in the DORA score function, as Steve is currently
# doing in his DORA optimization code.

# At this point, planning to put the DORA score in $tddr_dir/dora-score.
($tddr_dir, $responses_file) = @ARGV;

# Read the responses file, which tells us the directories where the several
# response Hill parameters can be found.
open(RESPONSES, $responses_file) or
    die("dora_score.pl: could not open "
	. "responses file $responses_file.");
while(<RESPONSES>)
{
    # In order to do any of these DORA scores right, it looks like
    # I will need to add the "nominal" population of each response
    # to the responses file.

    # This is needed to rescale the
    # baseline and range in the DORA score.  This will create havoc
    # in any situation where we're changing populations.  It means
    # connecting "responses" with actual species and getting those
    # species' "nominal" populations, which may be changed or not.
    #
    # An easier alternative, which makes nearly as much sense, is to
    # rescale differences between populations by the sum of those
    # populations.  A sort of "relative difference".
    chomp;
    ($dmp_file, $dmp_hdr) = split;
    $title = "$dmp_file--$dmp_hdr";
    push @titles, $title;
}
close(RESPONSES);

$totEC50 = 0.0;
$maxEC50 = POSIX::DBL_MIN;
$minEC50 = POSIX::DBL_MAX;

$totR0 = 0.0;
$maxR0 = POSIX::DBL_MIN;
$minR0 = POSIX::DBL_MAX;

$totAmp = 0.0;
$maxAmp = POSIX::DBL_MIN;
$minAmp = POSIX::DBL_MAX;

$maxN = POSIX::DBL_MIN;
$minN = POSIX::DBL_MAX;

# Mainly to kepp the number of "responses" as the length of the titles
# vector.
@rest_titles = @titles;
while(@rest_titles)
{
    # Construct the name of the Hill params file for this "response".
    $title = shift(@rest_titles);
    $hill_params_file = "$tddr_dir/$title.histos/hill-params";

    # Read the Hill parameters, accumulating them.
    open(PARAMS, $hill_params_file) or
	die("dora-score.pl: could not open params file $hill_params_file.");
    $_ = <PARAMS>;
    ($EC50, $R0, $RInf, $N) = split;
    close(PARAMS);

    $totEC50 += $EC50;
    if($maxEC50 < $EC50)
    {
	$maxEC50 = $EC50;
    }
    if($minEC50 > $EC50)
    {
	$minEC50 = $EC50;
    }

    $totR0 += $R0;
    if($maxR0 < $R0)
    {
	$maxR0 = $R0;
    }
    if($minR0 > $R0)
    {
	$minR0 = $R0;
    }

    $amp = $RInf - $R0;

    $totAmp += $amp;
    if($maxAmp < $amp)
    {
	$maxAmp = $amp;
    }
    if($minAmp > $amp)
    {
	$minAmp = $amp;
    }

    if($maxN < $N)
    {
	$maxN = $N;
    }
    if($minN > $N)
    {
	$minN = $N;
    }
}

# Weighting factors for the components of the dora score.  Subject to
# scrutiny and change.
$sigma_EC50 = 1.00;
$sigma_R0 = 1.00;
$sigma_Amp = 1.00;
$sigma_N = 1.00;

$EC50_spread = ($maxEC50 - $minEC50) / ($sigma_EC50 * (1.0 + $totEC50));
$R0_spread = ($maxR0 - $minR0) / ($sigma_R0 * (1.0 + $totR0));
$amp_spread = ($maxAmp - $minAmp) / ($sigma_Amp * (1.0 + $totAmp));
$N_spread = ($maxN - $minN) / $sigma_N;

$score_squared
    = ($EC50_spread * $EC50_spread)
    + ($R0_spread * $R0_spread)
    + ($amp_spread * $amp_spread)
    + ($N_spread * $N_spread);

$score = sqrt($score_squared);

# Write the DORA score to the output file.
$dora_score_file = "$tddr_dir/dora-score";
open(SCORE, ">$dora_score_file") or
    die("dora-score.pl: could not open output file $dora_score_file.");
print SCORE "$score\n";
close(SCORE);
