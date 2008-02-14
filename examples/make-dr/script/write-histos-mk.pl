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

($doses_file) = @ARGV;

# Length of the "tail" of the output to be taken as a sample of the
# equilibrium distriubiton, according to the MCMC principle.
$tail_lines = 200;

# The width of the histogram bins.
$histo_width = 5;

open(DOSES, $doses_file) or
    die("write-histos-mk.pl: could not open doses file $doses_file.");

# Skip the header line in the doses file.
<DOSES>;

while(<DOSES>)
{
    ($dose_dir, $dose) = split;

    print "\$(TDDR_HISTO_DIR)/$dose.cut : \$(TDDR_FILE)\n\n";

    print "\$(TDDR_HISTO_DIR)/$dose.data : \$(TDDR_HISTO_DIR)/$dose.cut\n";
    print "\ttail -n $tail_lines \$< > \$@\n\n";

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

print "\$(TDDR_HISTO_DIR)/all-drs : \$(TDDR_HISTO_DIR)/all-stats\n";
print "\t\$(SCRIPT)/normalize-dr.pl \$< \$@\n\n";

print "\$(TDDR_HISTO_DIR)/bezier : DR := \$(TDDR_HISTO_DIR)/bezier-dr\n";
print "\$(TDDR_HISTO_DIR)/bezier : DERIV := \$(TDDR_HISTO_DIR)/bezier-deriv\n";
print "\$(TDDR_HISTO_DIR)/bezier : VAR := \$(TDDR_HISTO_DIR)/bezier-var\n";
print "\$(TDDR_HISTO_DIR)/bezier :  \$(TDDR_HISTO_DIR)/all-stats\n";
print "\tnu -load \$(SCRIPT)/dose-distro.scm -stats-file \$< -dr-file \$(DR) -var-file \$(VAR) -deriv-file \$(DERIV)\n\n";

print "\$(TDDR_HISTO_DIR)/hill-params : \$(TDDR_HISTO_DIR)/all-stats\n";
print "\thillEst -sf \$< -of \$@\n\n";

print "\$(TDDR_HISTO_DIR)/histo-target : \\\n";
print "\t\$(TDDR_HISTO_DIR)/all-stats \\\n";
print "\t\$(TDDR_HISTO_DIR)/bezier \\\n";
print "\t\$(TDDR_HISTO_DIR)/hill-params \\\n";
print "\t\$(TDDR_HISTO_DIR)/all-drs\n\n";

close(DOSES);
