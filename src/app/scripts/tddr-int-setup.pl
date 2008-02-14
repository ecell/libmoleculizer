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

($input_file_path, $value_file, $sim_dir_stem) = @ARGV;

# Assuming that everything goes into the same directory
# as the value_file came from. 
($input_file_name, $input_file_dir, $suffix) = fileparse($input_file_path);

# Read line from the values file that tells how to do the substitution.
open(VALUE_HANDLE, $value_file);
$_ = <VALUE_HANDLE>;
@value_tokens = split;
close(VALUE_HANDLE);

# Generate list of values by parsing the expression.
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
	push @values, int($value);
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

open(BATCH_SCRIPT, ">$input_file_dir/batch_script.sh");
open(IMMED_SCRIPT, ">$input_file_dir/immed_script.sh");

# Using the value file's directory as the directory in which
# to generate simulation directories.  This parallels the report script.
($value_file_stem, $value_file_dir, $value_file_extension)
    = fileparse($value_file);
for($sim_ndx = 0;
    $sim_ndx < @values;
    $sim_ndx++)
{
    $sim_dir = "$value_file_dir/$sim_dir_stem\_$sim_ndx";
    mkdir $sim_dir;

    # Here is where tddr-float-setup.pl and tddr-int-setup.pl differ: in
    # tddr-float-setup.pl we use %g to write the value in floating point
    # notation, suitable for a reaction rate.  In tddr-int-setup.pl we use %d
    # to write the value as a decimal integer, suitable for a population.
    $value = sprintf "%d", ($values[$sim_ndx]);
    open(VALUE_HANDLE, ">$sim_dir/value");
    print VALUE_HANDLE "$value\n";
    close(VALUE_HANDLE);

    $varied_input = "$sim_dir/$input_file_name";

    system("java", "org.apache.xalan.xslt.Process",
	   "-in", "$input_file_path",
	   "-xsl", "$ENV{MOLECULIZER_DIR}/xml/xsl/edit-value.xsl",
	   "-param", "variable-marker", "VAR",
	   "-param", "new-value", "$value",
	   "-xml",
	   "-out", "$varied_input");

    print BATCH_SCRIPT "cd $sim_dir\n";
    print BATCH_SCRIPT "echo 'moleculizer < $input_file_name' | batch\n";
    print BATCH_SCRIPT "cd ..\n";

    print IMMED_SCRIPT "cd $sim_dir\n";
    print IMMED_SCRIPT "moleculizer < $input_file_name\n";
    print IMMED_SCRIPT "cd ..\n";

}

close(BATCH_SCRIPT);
close(IMMED_SCRIPT);

