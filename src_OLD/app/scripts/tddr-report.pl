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

($dmp_file, $value_file, $sim_dir_stem) = @ARGV;

# Parse the value_file to match up simulation directories with
# values.  This code is just copied from the setup script.
open(VALUE_HANDLE, $value_file);
$_ = <VALUE_HANDLE>;
@value_tokens = split;
close(VALUE_HANDLE);

# Parse the expression from the values file.
@values = ();
$operator = shift @value_tokens;
if($operator eq "+")
{
    ($value, $delta, $count) = @value_tokens;

    for($value_ndx = 0;
	$value_ndx < $count;
	$value_ndx++)
    {
	push @values, $value;
	$value = $value + $delta;
    }
}
elsif($operator eq "*")
{
    ($value, $factor, $count) = @value_tokens;

    for($value_ndx = 0;
	$value_ndx < $count;
	$value_ndx++)
    {
	push @values, $value;
	$value = $value * $factor;
    }
}
elsif($operator eq "=")
{
    @values = @value_tokens;
}
else
{
    die("unexpected operator $operator in $value_file.\n");
}

# Use the same recipe as the setup script for generating the
# simulation directory names.
($value_file_stem, $value_file_dir, $value_file_extension)
    = fileparse($value_file);

@sim_dir_list = ();
for($sim_ndx = 0;
    $sim_ndx < @values;
    $sim_ndx++)
{
    $sim_dir = "$value_file_dir/$sim_dir_stem\_$sim_ndx";
    push @sim_dir_list, $sim_dir;
}

# Read the headers out of the first .dmp file.
$typical_dmp_file = "$sim_dir_list[0]/$dmp_file";
open(DMP_FILE, $typical_dmp_file);
$_=<DMP_FILE>;
@col_headers = split;
close(DMP_FILE);

# Tease out the .dmp file name "stem," minus the .dmp,
# to construct the name of the tddr file.
($dmp_file_stem, $dmp_file_dir, $dmp_file_extension) 
    = fileparse($dmp_file, '\.dmp');

# Make a directory to hold the tddr reports for this file.
$tddr_dir = "tddrs";
mkdir $tddr_dir;

# Skip the time column, column 0.
for($col_ndx = 1;
    $col_ndx <= $#col_headers;
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
    while(@sim_dirs)
    {
	$sim_dir = shift @sim_dirs;

	# Cut the column of the .dmp file from this simulation.
	$tmp_cut_file = "/tmp/$$";
	$cut_ndx = $col_ndx + 1;
	system("cut -f $cut_ndx $sim_dir/$dmp_file > $tmp_cut_file");

	# Read the substituted value from the "value" file in the
	# simulation directory.
	open(VALUE_HANDLE, "$sim_dir/value");
	$substituted_value = <VALUE_HANDLE>;
	close(VALUE_HANDLE);
	chop $substituted_value;

	# Push the cut file onto the list of cut files
	# to be used in preparing the tddr file.
	$cut_file = "$cut_dir/$sim_dir";
	push @cut_files, ($cut_file);

	# Use the substituted value as a header in the cut file.
	open(CUT_FILE, ">$cut_file");
	print(CUT_FILE "$substituted_value\n");

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
    system("paste @cut_files > $tddr_file\n");

    # Remove the directory full of cuts.
    system("rm -r $cut_dir");
}
