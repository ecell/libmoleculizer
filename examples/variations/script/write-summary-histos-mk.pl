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

($doses_file, $variations_file, $output_dir) = @ARGV;

# The width of the histogram bins.
$histo_width = 5;

open(DOSES, $doses_file) or
    die("write-summary-histos-mk.pl: could not open doses file $doses_file.");

# Skip the header line in the doses file.
<DOSES>;

while(<DOSES>)
{
    ($dose_dir, $dose) = split;

    # We make the data file by cat'ing together all the data files from
    # the variation simulations.  Otherwise, everything should be more
    # or less the same.
    print "\$(TDDR_HISTO_DIR)/$dose.data : ";
    open(VARIATIONS, $variations_file) or
	die("write-summary-histos-mk.pl: could not open "
	    . "variations file $variations_file.");
    # Skip the header line in the variations file.
    <VARIATIONS>;
    while(<VARIATIONS>)
    {
	($variation_dir) = split;
	print "\\\n\t$output_dir/$variation_dir/tddrs/\$(DMP_FILE_STEM)--\$(COLUMN_HEADER).histos/$dose.data";
    }
    close(VARIATIONS);
    print "\n\tcat \$^ > \$@\n\n";

    print "\$(TDDR_HISTO_DIR)/$dose.histo : \$(TDDR_HISTO_DIR)/$dose.data\n";
    print "\trealbin $histo_width < \$< > \$@\n\n";

    print "\$(TDDR_HISTO_DIR)/$dose.stats : \$(TDDR_HISTO_DIR)/$dose.data\n";
    print "\t\$(SCRIPT)/stats.pl < \$< > \$@\n\n";

    print "\$(TDDR_HISTO_DIR)/$dose.gamma : \$(TDDR_HISTO_DIR)/$dose.stats\n";
    print "\tnu -load \$(SCRIPT)/gamma-density.scm < \$< > \$@\n\n";

    print "\$(TDDR_HISTO_DIR)/histo-target : \$(TDDR_HISTO_DIR)/$dose.histo \$(TDDR_HISTO_DIR)/$dose.gamma\n\n";
}

# Note that there is no real reason to do this here, it could just as well
# be in response.mk.
print "\$(TDDR_HISTO_DIR)/all-stats : \$(TDDR_HISTO_DIR)\n";
# Note that the at-sign needs to be escaped, since it is meaningful
# in Perl.
print "\t\$(SCRIPT)/write-all-stats.pl \$(\@D) $doses_file "
    . "\$(\@D)/all-stats\n\n";

print "\$(TDDR_HISTO_DIR)/histo-target : \$(TDDR_HISTO_DIR)/all-stats\n\n";

close(DOSES);
