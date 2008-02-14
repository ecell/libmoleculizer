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

# Take a "tail" of a Moleculizer .dmp output file, bin the
# columns into separte files, all with the same bin width, and use Gnuplot to
# display histograms of all the columns in a single plot.
#
# Also, one will want to extract the column headers of the .dmp file,
# so that they can be used to make the filenames for the binned data
# as well as the legends in the plot.

($input_file, $tail_lines, $bin_width, $output_directory) = @ARGV;

# Open the input file, and extract its first line as a vector, and close it
# again.
open(DATA_HANDLE, $input_file)
    or die("Could not open input file $input_file.");
$_ = <DATA_HANDLE>;
@column_headers = split;
close(DATA_HANDLE);

# Skip the first column, the simulation time with its undesirably funky
# column header.
shift @column_headers;

# For each of the remaining columns, generate a histogram file and
# arrange to plot the histogram file in the gnuplot script.
for($column_ndx = 0;
    $column_ndx < @column_headers;
    $column_ndx++)
{
    # Run cut to extract the corresponding column from the input file.
    $cut_field_ndx = $column_ndx + 2;
    $column_file = "$output_directory/$column_headers[$column_ndx].data";
    system("tail -n $tail_lines $input_file | cut -f$cut_field_ndx > $column_file");
}
