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

($dmp_file_stem, $column_header, $doses_file, $output_dir, $tddr_dir) = @ARGV;

# Open the doses file, so that we can use its first column, the simulation
# directories.
open(DOSES, $doses_file) or
    die("cut_dmp_column.pl: could not open doses file $doses_file.");

# Omit the header line of the doses file.
<DOSES>;

# Create a directory to hold cuts from the dump files.
$tmp_cut_dir = "$tddr_dir/tmp_cuts";
mkdir($tmp_cut_dir) or
    die("cut_dmp_column.pl: could not create directory $tmp_cut_dir.");

while(<DOSES>)
{
    ($sim_dir, $dose) = split;

    # The path to the dmp file produced by the simulation done in $sim_dir.
    $dmp_file_path = "$output_dir/$sim_dir/$dmp_file_stem.dmp";

    # The first time through, cut the first (sim-time) column out for
    # the tddr file.
    if(0 == @cut_files)
    {
	$cut_file = "$tmp_cut_dir/sim-time.cut";
	system("cut -f1 $dmp_file_path > $cut_file");
	push @cut_files, $cut_file;
    }

    # Open the dump file and read its header line.
    open(DMP, $dmp_file_path) or
	die("cut-dmp-column.pl: could not open dump file $dmp_file_path.");
    $_ = <DMP>;
    @header_tokens = split;
    close(DMP);

    # Laboriously find the desired field.
    $found_field_index = -1;
    $field_index = 1;
    while(@header_tokens)
    {
	$header_token = shift @header_tokens;
	if("$header_token" eq "$column_header")
	{
	    $found_field_index = $field_index;
	}
	$field_index++;
    }
    if($found_field_index < 0)
    {
	die("cut-dmp-column.pl: did not find column header "
	    . "$column_header in file $dmp_file_path.");
    }

    # Cut the desired field out of the .dmp file.
    $tmp_cut_file = "/tmp/cut-dmp-column-$$";
    system("cut -f$found_field_index $dmp_file_path > $tmp_cut_file");

    # Add the final cut file to the array of cut files, which we use
    # to paste all the cut files together in order.
    # $cut_file = "$tmp_cut_dir/$sim_dir.cut";
    $cut_file = "$tddr_dir/$dmp_file_stem--$column_header.histos/$dose.cut";
    push @cut_files, $cut_file;

    # Use the dose value as the header in the cut file.
    open(CUT_FILE, ">$cut_file") or
	die("cut-dmp-column.pl: could not open cut file $cut_file.");
    print CUT_FILE "$dose\n";

    # Copy everything except the header from the temporary cut file to
    # the final cut file.
    open(TMP_CUT_FILE, $tmp_cut_file) or
	die("cut-dmp-column.pl: could not open temporary "
	    . "cut file $tmp_cut_file.");
    # Skip the header.
    <TMP_CUT_FILE>;
    while(<TMP_CUT_FILE>)
    {
	print(CUT_FILE $_);
    }
    close(TMP_CUT_FILE);
    unlink $tmp_cut_file;
}
close(DOSES);

# Paste the cuts into the tddr file.
$tddr_file = "$tddr_dir/$dmp_file_stem--$column_header.tddr";
system("paste @cut_files > $tddr_file\n");

# Get rid of the temporary cuts.
system("rm -r $tmp_cut_dir");


