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

# This just replicates write-gp-histos-script.pl, but directs gnuplot to
# generate fig output instead of png.

# Writes gnuplot script to plot response histograms and fitted
# gamma densities.

($doses_file, $gp_script_file, $data_directory) = @ARGV;

print ("write-gp-doses-script: doses file is $doses_file.\n");
print ("write-gp-doses-script: script file is $gp_script_file.\n");
print ("write-gp-doses-script: data directory is $data_directory.\n");

# Read the doses, which appear as the second column in the doses
# file, as a vector.
@doses = ();
open(DOSES_HANDLE, $doses_file)
    or die("write-gp-doses-script.pl: could not open doses file $doses_file.");
# Skip the headers.
<DOSES_HANDLE>;
while(<DOSES_HANDLE>)
{
    ($dir_name, $dose) = split;
    push @doses, $dose;
}
close(DOSES_HANDLE);

# Write the gnuplot script.
open(GNUPLOT_HANDLE, ">$gp_script_file")
    or die("write-gp-doses-script.pl: could not open script file $gp_script_file.");

# Write the prelude of the gnuplot script, which does not depend on
# the column headers, which we want to appear in the legend.
print GNUPLOT_HANDLE "set term fig color\n";
print GNUPLOT_HANDLE "set yrange [0:*]\n";
print GNUPLOT_HANDLE "set data style histeps\n";
print GNUPLOT_HANDLE "set ylabel \"Fraction of sample points\"\n";
print GNUPLOT_HANDLE "set xlabel \"Number of molecules\"\n";
print GNUPLOT_HANDLE "plot ";

# Write plot items for all the histogram files and the gamma density files.
@rest_doses = @doses;
while(@rest_doses)
{
    $dose = shift @rest_doses;
    $histo_file = "$data_directory/$dose.histo";
    $gamma_density_file = "$data_directory/$dose.gamma";

    print GNUPLOT_HANDLE 
	"\"$histo_file\" title \"$dose\", ";
    print GNUPLOT_HANDLE 
 	"\"$gamma_density_file\" notitle lt 6";
   
    if(0 < @rest_doses)
    {
	print GNUPLOT_HANDLE ", ";
    }
    else
    {
	print GNUPLOT_HANDLE "\n";
    }
}
close(GNUPLOT_HANDLE);
