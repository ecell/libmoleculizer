#!/usr/bin/perl

$pi = 3.14159;
$liters_per_cubic_meter = 1.0e3;
$microns_per_meter = 1.0e6;

# Diameter of yeast, as of 29Jul04.
$diameter_in_microns = 5.0;

$diameter_in_meters = $diameter_in_microns / $microns_per_meter;

$radius_in_meters = $diameter_in_meters / 2.0;

# Assuming a yeast is spherical.
$volume_in_cubic_meters = (4.0 / 3.0)
    * $pi
    * $radius_in_meters
    * $radius_in_meters
    * $radius_in_meters;

$area_in_square_meters = 4.0 
    * $pi
    * $radius_in_meters
    * $radius_in_meters;

$volume_in_liters = $volume_in_cubic_meters * $liters_per_cubic_meter;

print "\nYeast shape:\n";
print "1) It has diameter about $diameter_in_microns microns.\n";
print "   (= $diameter_in_meters meters)\n";
print "2) It has a volume of about $volume_in_liters liters.\n";
print "3) It has surface area about $area_in_square_meters square meters.\n";

# Some of Kirsten's early populations.
$gpa1_count = 2390;
$ste4_count = 2045;
$ste20_count = 4212;
$ste5_count = 484;
$ste11_count = 3524;
$ste7_count = 924;
$fus3_count = 20413;
$kss1_count = 20795;
$ste12_count = 6016;
$far1_count = 1661;
$cdc42_count = 1595;
$pbs2_count = 2459;
$hog1_count = 5998;

print "\nKirsten's early populations:\n";
print "Gpa1\t";
print $gpa1_count;
print "\n";
print "Ste4\t";
print $ste4_count;
print "\n";
print "Ste20\t";
print $ste20_count;
print "\n";
print "Ste5\t";
print $ste5_count;
print "\n";
print "Ste11\t";
print $ste11_count;
print "\n";
print "Ste7\t";
print $ste7_count;
print "\n";
print "Fus3\t";
print $fus3_count;
print "\n";
print "Kss1\t";
print $kss1_count;
print "\n";
print "Ste12\t";
print $ste12_count;
print "\n";
print "Far1\t";
print $far1_count;
print "\n";
print "Cdc42\t";
print $cdc42_count;
print "\n";
print "Pbs2\t";
print $pbs2_count;
print "\n";
print "Hog1\t";
print $hog1_count;
print "\n";

print "\nYeast protein masses from MIPS database.\n";
print "(http://mips.gsf.de/projects/fungi)\n";
print "alpha\t";
print 13271.1;
print "\tdaltons\n";
print "No URL--look for an ad...\n\n";

print "Gpa1\t";
print 54075.9;
print "\tdaltons\n";
print "http://mips.gsf.de/genre/proj/yeast/searchEntryAction.do?text=YHR005c&db=\n\n";

print "Ste2\t";
print 47849.8;
print "\tdaltons\n";
print "http://mips.gsf.de/genre/proj/yeast/searchEntryAction.do?text=YFL026w&db=\n\n";

print "Ste4\t";
print 46581.2;
print "\tdaltons\n";
print "http://mips.gsf.de/genre/proj/yeast/searchEntryAction.do?text=YOR212w&db=\n\n";

print "Ste20\t";
print 102362;
print "\tdaltons\n";
print "http://mips.gsf.de/genre/proj/yeast/searchEntryAction.do?text=YHL007c&db=\n\n";

print "Ste5\t";
print 102728;
print "\tdaltons\n";
print "http://mips.gsf.de/genre/proj/yeast/searchEntryAction.do?text=YDR103w&db=\n\n";

print "Ste11\t";
print 80721.5;
print "\tdaltons\n";
print "http://mips.gsf.de/genre/proj/yeast/searchEntryAction.do?text=YLR362w&db=\n\n";

print "Ste7\t";
print 57709.4;
print "\tdaltons\n";
print "http://mips.gsf.de/genre/proj/yeast/searchEntryAction.do?text=YDL159w&db=\n\n";

print "Fus3\t";
print 40772.3;
print "\tdaltons\n";
print "http://mips.gsf.de/genre/proj/yeast/searchEntryAction.do?text=YBL016w&db=\n\n";

print "Kss1\t";
print 42692.5;
print "\tdaltons\n";
print "http://mips.gsf.de/genre/proj/yeast/searchEntryAction.do?text=YGR040w&db=\n\n";

print "Ste12\t";
print 77866.8;
print "\tdaltons\n";
print "http://mips.gsf.de/genre/proj/yeast/searchEntryAction.do?text=YHR084w&db=\n\n";

print "Far1\t";
print 94572.7;
print "\tdaltons\n";
print "http://mips.gsf.de/genre/proj/yeast/searchEntryAction.do?text=YJL157c&db=\n\n";

print "Cdc42\t";
print 21322.8;
print "\tdaltons\n";
print "http://mips.gsf.de/genre/proj/yeast/searchEntryAction.do?text=YLR229c&db=\n\n";

print "Pbs2\t";
print 72720.1;
print "\tdaltons\n";
print "http://mips.gsf.de/genre/proj/yeast/searchEntryAction.do?text=YJL128c&db=\n\n";

print "Hog1\t";
print 48858.6;
print "\tdaltons\n";
print "http://mips.gsf.de/genre/proj/yeast/searchEntryAction.do?text=YLR113w&db=\n\n";

print "\nOther useful masses.\n";
print "(From Wikipedia.)\n";
print "ATP\t";
print 507.181;
print "\tdaltons\n";
print "http://en.wikipedia.org/wiki/Adenosine_triphosphate\n\n";

print "ADP\t";
print 427.20;
print "\tdaltons\n";
print "http://en.wikipedia.org/wiki/Adenosine_diphosphate\n\n";

print "GTP\t";
print 523.18;
print "\tdaltons\n";
print "http://en.wikipedia.org/wiki/Guanosine_triphosphate\n\n";

print "GDP\t";
print 443.20;
print "\tdaltons\n";
print "http://en.wikipedia.org/wiki/Guanosine_diphosphate\n\n";

print "HPO_4\t";
print 95.97;
print "\tdaltons\n";
print "http://en.wikipedia.org/wiki/Phosphate#Chemical_properties\n\n";

# Typical diffusion rate of some E. coli proteins in square microns per
# second.
$elowitz_diffusion_number = 3.0;

$edn = $elowitz_diffusion_number
    / ($microns_per_meter * $microns_per_meter);

print "\nProteins in E. coli tend to diffuse at rates like\n";
print "$edn square meters per second.\n";
print "See Elowitz et al. 1999.\n\n";

print "\nKirsten says most on-rates of proteins are about 1.0e6. 19Oct06.\n\n";
