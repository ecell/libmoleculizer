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

use File::Basename;

($dmp_file, $value_file) = @ARGV;

# My new convention pertinent to this script is the following setup of
# columns of the value_file:

# DIR DOSE VAR1 VAR2 ...

# That is, the first column gives the simulation directory, relative to the
# current directory, the second column must give the value of the variable to
# be treated as a dose for purposes of this dose/response analysis.  My
# supposition now is that, if multiple simulations use the same dose, then
# that dose will appear more than once in the reports.

# Parse the first two columns of the value file, which give the simulation
# directories and the doses.
open(VALUE_HANDLE, $value_file)
    or die("new-report.pl: could not open values file $value_file.");

# Skip the headers in the values file.
<VALUE_HANDLE>;

# Process the rest of the lines, which give directories and doses.
@sim_dir_list = ();
@doses_list = ();
while(<VALUE_HANDLE>)
{
    @value_line_tokens = split;
    
    $sim_dir = shift @value_line_tokens;
    push @sim_dir_list, $sim_dir;

    $dose = shift @value_line_tokens;
    push @doses_list, $dose;
}
close(VALUE_HANDLE);

# Read the headers out of the first .dmp file.
$typical_dmp_file = "$sim_dir_list[0]/$dmp_file";
open(DMP_FILE, $typical_dmp_file)
    or die("new-report.pl: could not open dmp file $typical_dmp_file.");
$_=<DMP_FILE>;
@col_headers = split;
close(DMP_FILE);

print "new-report.pl: column headers are ", @col_headers, "\n";

# Tease out the .dmp file name "stem," minus the .dmp,
# to construct the name of the tddr file.
($dmp_file_stem, $dmp_file_dir, $dmp_file_extension) 
    = fileparse($dmp_file, '\.dmp');

# Make a directory to hold the tddr reports for this file.
$tddr_dir = "tddrs";
mkdir $tddr_dir;

# Skip the time column, column 0.
for($col_ndx = 1;
    $col_ndx < @col_headers;
    $col_ndx++)
{
    # Make directory to hold the cuts.
    $cut_dir = "$tddr_dir/$dmp_file_stem--cuts";
    mkdir $cut_dir;

    # Generate the simulation time cut.
    system("cut -f 1 $sim_dir_list[0]/$dmp_file > $cut_dir/simulation_time");

    # Write cut file for each .dmp file.
    @cut_files = ("$cut_dir/simulation_time");
    @sim_dirs = @sim_dir_list;
    @doses = @doses_list;
    while(@sim_dirs)
    {
	$sim_dir = shift @sim_dirs;
	$dose = shift @doses;

	# Cut the column of the .dmp file from this simulation.
	$tmp_cut_file = "/tmp/$$";
	$cut_ndx = $col_ndx + 1;
	system("cut -f $cut_ndx $sim_dir/$dmp_file > $tmp_cut_file");

	# Push the cut file onto the list of cut files
	# to be used in preparing the tddr file.
	$cut_file = "$cut_dir/$sim_dir";
	push @cut_files, ($cut_file);

	# Use the dose value as a header in the cut file.
	open(CUT_FILE, ">$cut_file");
	print(CUT_FILE "$dose\n");

	# Copy everything except the header from the temporary cut file.
	open(TMP_CUT_FILE, "$tmp_cut_file");
	<TMP_CUT_FILE>;
	while(<TMP_CUT_FILE>)
	{
	    print(CUT_FILE $_);
	}
	close(TMP_CUT_FILE);
	unlink ($tmp_cut_file);

	close(CUT_FILE);
    }
    
    # Paste the cuts into the tddr file.
    $tddr_file = "$tddr_dir/$dmp_file_stem--$col_headers[$col_ndx].tddr";
    print "new-report.pl: writing tddr file $tddr_file.\n";
    system("paste @cut_files > $tddr_file\n");

    # Remove the directory full of cuts.
    system("rm -r $cut_dir");
}
