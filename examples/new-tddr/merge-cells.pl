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

use File::Copy;
use File::Basename;

# Parse the command line.
($variations_file_path, $doses_file_path) = @ARGV;

# Parse the file names.
($variations_file_name, $variations_file_dir, $variations_file_suffix)
    = fileparse($variations_file_path);
($doses_file_name, $doses_file_dir, $doses_file_suffix)
    = fileparse($doses_file_path);

# For each dose, merge the "tails" of the responses from each cell type
# into a single large file of the same type for that dose.

# Read the doses file, assembling the doses
# in a vector, so that we don't have to read this file over again for
# each cell type.  The doses themselves are used as file names, a dubious
# practice, which I may want to correct.
open(DOSES, $doses_file_path)
    or die("merge-cells.pl: could not open doses file $doses_file_path.");
# Skip the header line.
<DOSES>;
@doses = ();
while(<DOSES>)
{
    # Parse the row of the doses file into fields.
    @dose_entries = split;

    # The first column gives the directory for this dose, which we don't
    # need at this point.
    shift @dose_entries;

    # The second column gives the dose, which is actually used to
    # to make the file name of the sample for that dose.  Reusing the
    # directory name might have been a better choice.
    $dose = shift @dose_entries;
    push @doses, $dose;
}
close(DOSES);

# Some setup-based constants.
$output_dir = "new-tddr.out";

# We expect this output directory, analogous to those in the separate
# cell-type directories, to be created by the main makefile.
$target_histo_dir = "$output_dir/histos";

open(VARIATIONS, $variations_file_path)
    or die("merge-cells.pl: could not open variations file $variations_file_path.");
# Omit the header line.
<VARIATIONS>;
while(<VARIATIONS>)
{
    # This is the directory, followed by the parameter variations for this
    # cell type.
    @variation = split;

    # We only need the directory.
    $variation_dir = shift @variation;
    $variation_histo_dir = "$variation_dir/histos";

    # For each dose, merge the samples of this cell type (i.e. variation) with
    # the samples of the other cell types.
    @rest_doses = @doses;
    while(@rest_doses)
    {
	$dose = shift @rest_doses;

	# Cook up the names of the sample file for this cell type
	# and of the main sample file.
	$dose_file_name = "$dose.data";
	$type_dose_file = "$variation_histo_dir/$dose_file_name";
	$main_dose_file = "$target_histo_dir/$dose_file_name";

	print("merge-cells.pl: merging $type_dose_file ");
	print("into $main_dose_file.\n");

	# Merge the sample for this cell type into the
	# main sample.
	system("cat $type_dose_file >> $main_dose_file");
    }
}
close(VARIATIONS);
