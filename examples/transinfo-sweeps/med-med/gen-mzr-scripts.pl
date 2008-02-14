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

# This program generates scripts that run in a single dose/response simulation.

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

# The first line of the values file gives the strings (like VAR1, VAR2,
# etc.) that are to be substituted in the input file with values on remaining
# lines of the tab_values_file.  The first column in this special line is
# special, it's a dummy placeholder for which I conventionally use DIR.
open(VALUE_HANDLE, $values_file_path);
$_ = <VALUE_HANDLE>;
@sub_targets = split;

open(BATCH_SCRIPT, ">batch_script.sh");
open(IMMED_SCRIPT, ">immed_script.sh");

while(<VALUE_HANDLE>)
{
    # All the columns of a particular substitution.
    @sub_values = split;

    # Again, the first column, which gives a relative path to the simulation
    # output directory, is special.
    $sim_dir_relative_path = shift @sub_values;

    # The job for each line in the values file is executed in a directory
    # whose path relative to the VALUES FILE (not the input file) is given in
    # the first column of the line.
    $sim_dir = "$values_file_dir/$sim_dir_relative_path";

    # Make addition to the script that runs the simulations in batch mode.
    print BATCH_SCRIPT "cd $sim_dir\n";
    print BATCH_SCRIPT "echo 'moleculizer < $input_file_name' | batch\n";
    print BATCH_SCRIPT "cd ..\n";

    # Make addition to the script that runs the simulations now.
    print IMMED_SCRIPT "cd $sim_dir\n";
    print IMMED_SCRIPT "moleculizer < $input_file_name\n";
    print IMMED_SCRIPT "cd ..\n";
}

close(BATCH_SCRIPT);
close(IMMED_SCRIPT);


