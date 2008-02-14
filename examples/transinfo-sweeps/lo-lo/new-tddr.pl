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
($input_file_path, $values_file_path) = @ARGV;

# Parse the file names.  

# The job for each line in the values file is executed in a directory whose
# path relative to the VALUES FILE (not the input file) is given in the first
# column of the line.
($input_file_name, $input_file_dir, $input_file_suffix)
    = fileparse($input_file_path);
($values_file_name, $values_file_dir, $values_file_suffix)
    = fileparse($values_file_path);

# The first line of the tab_values_file gives the property values (like
# VAR1, VAR2, etc.) that are to be substituted with values on remaining lines
# of the tab_values_file.
open(VALUE_HANDLE, $values_file_path);
$_ = <VALUE_HANDLE>;
@sub_targets = split;

# The first header, on the output directory column, is a dummy, by my
# convention DIR.  The remaining strings in this array are to be replaced
# in the input file with the values on the columns below them.
shift @sub_targets;

# Just a part of doing multiple substitutions.
$substitution_source_path = "/tmp/"."$$"."_tddr_sub_source";
$substitution_target_path = "/tmp/"."$$"."_tddr_sub_target";

# Now read one line at a time from the tab_values_file to get values to
# substitute for the "symbolic" property values (like VAR1, VAR2, etc.)
# given on the first line of the file.
while(<VALUE_HANDLE>)
{
    # To begin the loop, copy the input file to the substitution source.
    copy($input_file_path, $substitution_source_path);

    @sub_values = split;

    # Again, the first column, which gives a relative path to the simulation
    # output directory, is special.
    $sim_dir_relative_path = shift @sub_values;

    # Perform the substitutions one at a time.
    for($sub_ndx = 0;
	$sub_ndx < @sub_targets;
	$sub_ndx++)
    {
	# Do one substitution.
	system("java", "org.apache.xalan.xslt.Process",
	       "-in", "$substitution_source_path",
	       "-xsl", "$ENV{MOLECULIZER_DIR}/xml/xsl/edit-value.xsl",
	       "-param", "variable-marker", "$sub_targets[$sub_ndx]",
	       "-param", "new-value", "$sub_values[$sub_ndx]",
	       "-xml",
	       "-out", "$substitution_target_path");

	# Move the substitution_target_path to the substitution_source_path
	# for the next iteration.
	unlink $substitution_source_path;
	copy($substitution_target_path, $substitution_source_path);
	unlink $substitution_target_path;
    }

    # Make the simulation output directory. The job for each line in the
    # values file is executed in a directory whose path relative to the VALUES
    # FILE (not the input file) is given in the first column of the line.
    $sim_dir = "$values_file_dir/$sim_dir_relative_path";
    print "new-tddr.pl creating directory $sim_dir\n";
    mkdir $sim_dir
	or die("Could not make simulation output directory $sim_dir.");

    # Copy the simulation input file, with substitutions made, to the
    # simulation directory.
    copy($substitution_source_path, "$sim_dir/$input_file_name");
    unlink $substitution_source_path;
}
