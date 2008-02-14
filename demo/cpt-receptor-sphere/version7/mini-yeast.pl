#!/usr/bin/perl

# Factor by which linear size is reduced.
$mini_factor = 0.5;

# Factor for reducing the number of molecules of different species.
$mini_vol_factor = $mini_factor * $mini_factor * $mini_factor;

$pi = 3.14159;
$liters_per_cubic_meter = 1.0e3;
$microns_per_meter = 1.0e6;

# Diameter of yeast, as of 29Jul04.
$diameter_in_microns = 5.0 * $mini_factor;

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

print "\nMini-yeast is $mini_factor times the (linear) size\n";
print "of a regular yeast.\n";

print "\nMini-yeast shape:\n";
print "1) It has diameter about $diameter_in_microns microns.\n";
print "   (= $diameter_in_meters meters.)\n";
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

print "\nMini-yeast populations derived from\n";
print "Kirsten's early numbers:\n";
print "Gpa1\t";
print $gpa1_count * $mini_vol_factor;
print "\n";

print "Ste4\t";
print $ste4_count * $mini_vol_factor;
print "\n";

print "Ste20\t";
print $ste20_count * $mini_vol_factor;
print "\n";

print "Ste5\t";
print $ste5_count * $mini_vol_factor;
print "\n";

print "Ste11\t";
print $ste11_count * $mini_vol_factor;
print "\n";

print "Ste7\t";
print $ste7_count * $mini_vol_factor;
print "\n";

print "Fus3\t";
print $fus3_count * $mini_vol_factor;
print "\n";

print "Kss1\t";
print $kss1_count * $mini_vol_factor;
print "\n";

print "Ste12\t";
print $ste12_count * $mini_vol_factor;
print "\n";

print "Far1\t";
print $far1_count * $mini_vol_factor;
print "\n";

print "Cdc42\t";
print $cdc42_count * $mini_vol_factor;
print "\n";

print "Pbs2\t";
print $pbs2_count * $mini_vol_factor;
print "\n";

print "Hog1\t";
print $hog1_count * $mini_vol_factor;
print "\n";

print "\nYeast protein masses from MIPS database.\n";
print "(http://mips.gsf.de/projects/fungi)\n";
print "alpha\t";
print 13271.1;
print "\tdaltons\n";

print "Gpa1\t";
print 54075.9;
print "\tdaltons\n";

print "Ste2\t";
print 47849.8;
print "\tdaltons\n";

print "Ste4\t";
print 46581.2;
print "\tdaltons\n";

print "Ste20\t";
print 102362;
print "\tdaltons\n";

print "Ste5\t";
print 102728;
print "\tdaltons\n";

print "Ste11\t";
print 80721.5;
print "\tdaltons\n";

print "Ste7\t";
print 57709.4;
print "\tdaltons\n";

print "Fus3\t";
print 40772.3;
print "\tdaltons\n";

print "Kss1\t";
print 42692.5;
print "\tdaltons\n";

print "Ste12\t";
print 77866.8;
print "\tdaltons\n";

print "Far1\t";
print 94572.7;
print "\tdaltons\n";

print "Cdc42\t";
print 21322.8;
print "\tdaltons\n";

print "Pbs2\t";
print 72720.1;
print "\tdaltons\n";

print "Hog1\t";
print 48858.6;
print "\tdaltons\n";

print "\nOther useful masses.\n";
print "(From Wikipedia.)\n";
print "ATP\t";
print 507.181;
print "\tdaltons\n";

print "ADP\t";
print 427.20;
print "\tdaltons\n";

print "GTP\t";
print 523.18;
print "\tdaltons\n";

print "GDP\t";
print 443.20;
print "\tdaltons\n";

print "HPO_4\t";
print 95.97;
print "\tdaltons\n";

$cytoplasm_compartment_count = 17;
$membrane_compartment_count = 16;

# The ratio of a membrane compartment volume to a cytoplasm compartment volume.
$membrane_size_factor = 0.1;

# Multiply cell volume by this to get cytoplasm volume.
# Also useful for getting cytoplasm compartment populations.
$cytoplasm_volume_factor = 1.0 / ($cytoplasm_compartment_count
				  + ($membrane_size_factor *
				     $membrane_compartment_count));

# Multiply cell volume by this to get membrane volume.
# Also useful for getting membrane compartment populations.
$membrane_volume_factor = $membrane_size_factor * $cytoplasm_volume_factor;

$cytoplasm_compartment_volume
    = $volume_in_liters * $cytoplasm_volume_factor;

$membrane_compartment_volume
    = $volume_in_liters * $membrane_volume_factor;

# Area exposed to the outside and to membrane compartments.
$membrane_exterior_area
    = $area_in_square_meters / $membrane_compartment_count;

# Area between membrane compartments.
$membrane_adjacent_area
    = $membrane_exterior_area * $membrane_size_factor;

print "\nAssuming that there are\n";
print $cytoplasm_compartment_count;
print "\ncytoplasm compartments,and that each of the\n";
print $membrane_compartment_count;
print "\nmembrane compartments is\n";
print $membrane_size_factor;
print "\ntimes the volume of a cytoplasm compartment,\n";
print "then each cytoplasm compartment should be\n";
print $cytoplasm_compartment_volume;
print "\nliters, and each membrane compartment should be\n";
print $membrane_compartment_volume;
print "\nliters.\n";

print "\nEach membrane compartment should have exterior/interior area\n";
print $membrane_exterior_area;
print "\nsquare meters.\n";
print "Boundaries between membrane compartments should be\n";
print $membrane_adjacent_area;
print "\nsquare meters.\n";

print "\nTrial distance from membrane compartment to underlying cytoplasm\n";
print "compartment is\n";
print $radius_in_meters / 3.0;
print "\nmeters.  Trial distance between cytoplasm compartments\n";
print "and between membrane compartments is\n";
print 2.0 * $radius_in_meters / 3.0;
print "\nmeters.\n";

print "\nMini-cell compartment pops.\n";
print "Protein\tcytoplasm cpt\t\tmembrane cpt\n";
print "Gpa1\t";
print $gpa1_count * $mini_vol_factor * $cytoplasm_volume_factor;
print "\t";
print $gpa1_count * $mini_vol_factor * $membrane_volume_factor;
print "\n";

print "Ste4\t";
print $ste4_count * $mini_vol_factor * $cytoplasm_volume_factor;
print "\t";
print $ste4_count * $mini_vol_factor * $membrane_volume_factor;
print "\n";

print "Ste20\t";
print $ste20_count * $mini_vol_factor * $cytoplasm_volume_factor;
print "\t";
print $ste20_count * $mini_vol_factor * $membrane_volume_factor;
print "\n";

print "Ste5\t";
print $ste5_count * $mini_vol_factor * $cytoplasm_volume_factor;
print "\t";
print $ste5_count * $mini_vol_factor * $membrane_volume_factor;
print "\n";

print "Ste11\t";
print $ste11_count * $mini_vol_factor * $cytoplasm_volume_factor;
print "\t";
print $ste11_count * $mini_vol_factor * $membrane_volume_factor;
print "\n";

print "Ste7\t";
print $ste7_count * $mini_vol_factor * $cytoplasm_volume_factor;
print "\t";
print $ste7_count * $mini_vol_factor * $membrane_volume_factor;
print "\n";

print "Fus3\t";
print $fus3_count * $mini_vol_factor * $cytoplasm_volume_factor;
print "\t";
print $fus3_count * $mini_vol_factor * $membrane_volume_factor;
print "\n";

print "Kss1\t";
print $kss1_count * $mini_vol_factor * $cytoplasm_volume_factor;
print "\t";
print $kss1_count * $mini_vol_factor * $membrane_volume_factor;
print "\n";

print "Ste12\t";
print $ste12_count * $mini_vol_factor * $cytoplasm_volume_factor;
print "\t";
print $ste12_count * $mini_vol_factor * $membrane_volume_factor;
print "\n";

print "Far1\t";
print $far1_count * $mini_vol_factor * $cytoplasm_volume_factor;
print "\t";
print $far1_count * $mini_vol_factor * $membrane_volume_factor;
print "\n";

print "Cdc42\t";
print $cdc42_count * $mini_vol_factor * $cytoplasm_volume_factor;
print "\t";
print $cdc42_count * $mini_vol_factor * $membrane_volume_factor;
print "\n";

print "Pbs2\t";
print $pbs2_count * $mini_vol_factor * $cytoplasm_volume_factor;
print "\t";
print $pbs2_count * $mini_vol_factor * $membrane_volume_factor;
print "\n";

print "Hog1\t";
print $hog1_count * $mini_vol_factor * $cytoplasm_volume_factor;
print "\t";
print $hog1_count * $mini_vol_factor * $membrane_volume_factor;
print "\n";

# Typical diffusion rate of some E. coli CYTOPLASMIC proteins in square
# microns per second.  
$elowitz_diffusion_number = 3.0;

$edn = $elowitz_diffusion_number
    / ($microns_per_meter * $microns_per_meter);

print "\nCytoplasmic proteins in E. coli tend to diffuse at rates like\n";
print "$edn square meters per second.\n";
print "See ";

print "Elowitz, M., Surette, M., Wolf P-E.  Stock, J.,and Leibler, S.\n";
print "Protein mobility in the cytoplasm of Escherichia coli.\n";
print "Journal of Bacteriology, 181,197-203 (1999).\n";

# Diffusion rate of some E. coli MEMBRANE proteins in square microns per
# second.
$deich_moerner_ecoli = 12.0e-3;

$dme = $deich_moerner_ecoli
    / ($microns_per_meter * $microns_per_meter);

print "\nMembrane proteins in E. coli tend to diffuse at rates like\n";
print "$dme square meters per second.\n";
print "See ";

print "Deich, J., Judd, E., McAdams, H. and Moerner, W.\n";
print "Visualization of the movement of single histidine\n";
print "kinase molecules in live Caulobacter cells.\n";
print "PNAS 101,15921-15926 (2004).\n";

# Diffusion rate of some EUKARYOTIC MEMBRANE proteins in square microns per
# second.
$deich_moerner_euk_low = 5.0e-3;
$deich_moerner_euk_hi = 500.0e-3;

$dmk_lo = $deich_moerner_euk_low
    / ($microns_per_meter * $microns_per_meter);

$dmk_hi = $deich_moerner_euk_hi
    / ($microns_per_meter * $microns_per_meter);

print "\nMembrane proteins in eukaryotes tend to diffuse at rates\n";
print "in the range\n";
print "$dmk_lo - $dmk_hi\n";
print "square meters per second.\n";
print "See Deich et al. above.\n";

# Trial diffusion rate for small molecules (yes, all of them.)
# Yes, just a fraction of the rate for cytoplasmic proteins.
$small_molecule_factor = 1.0e3;
print "\nTrial diffusion rate for small molecules:\n";
print $small_molecule_factor * $edn;
print "\nsquare meters per second.\n";

print "\n";

# Factor by which the rates of binary reactions were reduced to compensate
# for reduced volume
print "Binary reaction rates reduced by factor of 10 to compensate for\n";
print "reduced volume.\n";
    
