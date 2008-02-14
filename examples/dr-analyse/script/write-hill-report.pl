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

($tddr_dir, $responses_file, $report_file) = @ARGV;

open(RESPONSES, $responses_file) or
    die("write-hill-report.pl: could not open "
	. "responses file $responses_file.");
while(<RESPONSES>)
{
    chomp;
    ($dmp_file, $dmp_hdr) = split;
    $title = "$dmp_file--$dmp_hdr";
    push @titles, $title;

}
close(RESPONSES);

open(REPORT, ">$report_file") or
    die("write-hill-report.pl: could not open "
	. "output file $report_file.");

@rest_titles = @titles;
while(@rest_titles)
{
    $title = shift @rest_titles;
    $hill_params_file = "$tddr_dir/$title.histos/hill-params";

    open(PARAMS, "$hill_params_file") or
	die("write-hill-report.pl: could not open "
	    . "hill params file $hill_params_file.");
    chomp($params = <PARAMS>);
    close(PARAMS);

    print REPORT "$title\t$params\n";
}
close(REPORT);
