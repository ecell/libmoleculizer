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

($tddr_file, $timepoint) = @ARGV;

# Open the tddr file, produced by tddr-report.pl.
open(TDDR_HANDLE, $tddr_file)
    or die("Could not open tddr file $tddr_file.\n");

# Read the header line from the tddr file. The tokens after the first
# on this line should give us the doses.
$_ = <TDDR_HANDLE>;
@header_tokens = split;

# Look for the timepoint, as a floating point number, as the first token on
# a following line.
while(<TDDR_HANDLE>)
{
    @response_tokens = split;
    last if($response_tokens[0] == $timepoint);
}
close(TDDR_HANDLE);

# Throw out the first header token, which is the header token for the
# simulation time.
shift @header_tokens;

# Throw out the first token from the target line, which is the simulation time.
$test_time = shift @response_tokens;

# Print the rest to a temporary file as dose/response pairs.
# We will plot this temporary file with gnuplot.
$tmp_dr_file = "/tmp/tddr-plot-timepoint-$$";
open(TMP_DR_HANDLE, ">$tmp_dr_file")
    or die("Could not open temporary file $tmp_dr_file for writing.\n");

# By reporting the time from the line on the plot, we will always check that
# we're getting the timepoint we expect.
while(@header_tokens)
{
    print(TMP_DR_HANDLE
	  shift(@header_tokens),
	  "\t",
	  shift(@response_tokens),
	  "\n");
}
close(TMP_DR_HANDLE);

# Print the name of the temporary file, so that we can plot it from a shell
# script.
print "$tmp_dr_file\n";
