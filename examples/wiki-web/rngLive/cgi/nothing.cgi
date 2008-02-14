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

# Test cgi document.

use CGI;

$working_dir = $ENV{'PWD'};

$path = $ENV{'PATH'};

$test_file_name = "../foo";

open(TEST_HANDLE, $test_file_name)
    or die("could not open test file $test_file_name.");

$query = new CGI;

print $query->header;

print $query->start_html(-title=>'Simple test');

print $query->h1("A simple test");

print $query->hr;

print $query->p("The working directory is $working_dir.");
print $query->p("PATH is $path.");
if($_ = <TEST_HANDLE>)
{
    print "The test file contained $_\n";
}
else
{
    print "Couldn't read from test file.\n";
}

print $query->hr;

print $query->end_html;




