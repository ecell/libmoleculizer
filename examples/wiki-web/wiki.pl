#!/usr/bin/perl

# This script enshrines the numbers from Ty's simulation, as conveyed to me
# by Kirsten, Wed Oct 25, and converts them into forms/units useable by
# Moleculizer.
print "\n";
print "Conversions of units and other computations to get Moleculizer\n";
print "parameters from Wiki parameters and parameters from Ty's";
print "simulation.\n\n";

# Units conversion factors.
$litersPerFemtoliter = 1.0e-15;
$molesPerMicromole = 1.0e-6;

# Some standard rates.
$diffusionLimitedOnRate = 1.0e6;
$noBindingOnRate = 1.0;
$noBindingOffRate = 1.0e10;
$noUnbindingOffRate = 1.0e-10;

print "-----------\n\n";

print "Nucleotide binding and unbinding by Gpa1.\n";
print "(Not involved in nuclotide exchange.)\n\n";

print "Gpa1 binding to GDP.\n";
$Gpa1GdpPreferenceFactor = 10;
$kon_Gpa1_GDP = $Gpa1GdpPreferenceFactor * $diffusionLimitedOnRate;
$koff_Gpa1_GDP = $noUnbindingOffRate;
print "kon_Gpa1_GDP = $kon_Gpa1_GDP HzPerM\n";
print "koff_Gpa1_GDP = $koff_Gpa1_GDP Hz\n\n";

print "Gpa1 binding to GTP.\n";
$kon_Gpa1_GTP = $diffusionLimitedOnRate;
$koff_Gpa1_GTP = $noUnbindingOffRate;
print "kon_Gpa1_GTP = $kon_Gpa1_GTP HzPerM\n";
print "koff_Gpa1_GTP = $koff_Gpa1_GTP Hz\n\n";

print "-----------\n\n";

print "Nucleotide exchange and hydrolysis.\n\n";

# Factor by which binding to GDP enhances binding of Gpa1 to Ste4,
# presumptively by reducing the off-rate. Description refers only to
# dissociation constant.
#
# This is also the factor by which the rate of exchange of GDP for GTP is
# decreased when Gpa1 is bound to Ste4, regardless of whether Gpa1 is bound to
# Ste2.
$GTP_modulating_factor = 1000;

print "Autocatalyzed exchange of GDP for GTP by Gpa1. (Gpa1 not bound to\n";
print "Ste4-Ste18.)\n";
$kcat_Gpa1GDP_GTP = 6.17e-4;
print "kcat_Gpa1GDP_GTP = $kcat_Gpa1GDP_GTP Hz\n\n";

print "Autocatalyzed exchange of GDP for GTP by Gpa1 bound to Ste4-Ste18.\n";
$kcat_Gpa1GDPSte4Ste18_GTP = $kcat_Gpa1GDP_GTP / $GTP_modulating_factor;
print "kcat_Gpa1GDPSte4Ste18_GTP = $kcat_Gpa1GDPSte4Ste18_GTP Hz\n\n";

print "alpha-Ste2 catalyzed exchange of GDP for GTP\n";
print "by Gpa1 bound to Ste4-Ste18.\n";
#
# Perhaps take this to be the off-rate of GDP when in this configuration,
# assuming plenty of GTP and a reasonably fast on-rate.
#
$kcat_PheromoneSte2Gpa1GDPSte4Ste18_GTP = 1.5;
print "kcat_PheromoneSte2Gpa1GDPSte4Ste18_GTP = $kcat_PheromoneSte2Gpa1GDPSte4Ste18_GTP Hz\n\n";

print "alpha-Ste2 catalyzed exchange of GDP for GTP by Gpa1. (Gpa1 \n";
print "not bound to Ste4.)\n";
#
# Perhaps take this to be the off-rate of GDP when in this configuration,
# assuming plenty of GTP and a reasonably fast on-rate.
#
$kcat_PheromoneSte2Gpa1GDP_GTP
    = $GTP_modulating_factor * $kcat_PheromoneSte2Gpa1GDPSte4Ste18_GTP;
print "kcat_PheromoneSte2Gpa1GDP_GTP = $kcat_PheromoneSte2Gpa1GDP_GTP Hz\n\n";

print "Autocatalyzed hydrolysis of GTP by Gpa1.\n";
$kcat_Gpa1GTP_GDP = 5.0e-3;
print "kcat_Gpa1GTP_GDP = $kcat_Gpa1GTP_GDP Hz\n\n";

print "Sst2 catalyzed hydrolysis of GTP-Gpa1.\n";
$kcat_Sst2Gpa1GTP_GDP = 4.0;
print "kcat_Sst2Gpa1GTP_GDP = $kcat_Sst2Gpa1GTP_GDP Hz.\n\n";

print "-----------\n\n";

print "Binding of alpha factor to Ste2.\n";
# As supplied.
$Kd_Pheromone_Ste2_micromoles = 4.9e-3;
$Kd_Pheromone_Ste2 = $Kd_Pheromone_Ste2_micromoles * $molesPerMicromole;
$koff_Pheromone_Ste2 = 1.3e-3;
$kon_Pheromone_Ste2 = $koff_Pheromone_Ste2 / $Kd_Pheromone_Ste2;
print "kon_Pheromone_Ste2 = $kon_Pheromone_Ste2\n";
print "koff_Pheromone_Ste2 = $koff_Pheromone_Ste2\n\n";

print "-----------\n\n";

print "Binding of Ste2 to Gpa1.\n\n";

print "Binding of Ste2 to Gpa1 that is not bound to Ste4.\n";
# As provided.
$kon_Ste2_Gpa1_HzPerMicromole = 20;
$kon_Ste2_Gpa1 = $kon_Ste2_Gpa1_HzPerMicromole / $molesPerMicromole;
$koff_Ste2_Gpa1 = 20;
print "kon_Ste2_Gpa1 = $kon_Ste2_Gpa1 HzPerM\n";
print "koff_Ste2_Gpa1 = $koff_Ste2_Gpa1 Hz\n\n";

# Multiplier by which binding to Ste4 increases Gpa1's affinity for Ste2.
# By doctrine, this should also be the factor by which binding to Ste2 should
# increase Gpa1's affinity for Ste4.
$Ste2_Gpa1_Ste4Ste18_coop_factor = 10;

print "Binding of Ste2 to Ste4-bound Gpa1.\n";
$kon_Ste2_Gpa1Ste4Ste18 = $kon_Ste2_Gpa1;
$koff_Ste2_Gpa1Ste4Ste18 = $koff_Ste2_Gpa1 / $Ste2_Gpa1_Ste4Ste18_coop_factor;
print "kon_Ste2_Gpa1Ste4Ste18 = $kon_Ste2_Gpa1Ste4Ste18 HzPerM\n";
print "koff_Ste2_Gpa1Ste4Ste18 = $koff_Ste2_Gpa1Ste4Ste18 Hz\n\n";

print "-----------\n\n";

print "Binding of Gpa1 with Ste4-18.\n\n";

print "Binding of Gpa1-GDP with Ste4-Ste18. (Gpa1 not bound to Ste2.)\n";
# As provided.
$kon_Gpa1GDP_Ste4Ste18_HzPerMicromole = 50;
$kon_Gpa1GDP_Ste4Ste18
    = $kon_Gpa1GDP_Ste4Ste18_HzPerMicromole / $molesPerMicromole;
$koff_Gpa1GDP_Ste4Ste18 = 0.5;
print "kon_Gpa1GDP_Ste4Ste18 = $kon_Gpa1GDP_Ste4Ste18 HzPerM\n";
print "koff_Gpa1GDP_Ste4Ste18 = $koff_Gpa1GDP_Ste4Ste18 Hz\n\n";

print "Binding of Gpa1-GTP with Ste4-Ste18. (Gpa1 not bound to Ste2.)\n";
$kon_Gpa1GTP_Ste4Ste18 = $kon_Gpa1GDP_Ste4Ste18;
$koff_Gpa1GTP_Ste4Ste18 = $GTP_modulating_factor * $koff_Gpa1GDP_Ste4Ste18;
print "kon_Gpa1GTP_Ste4Ste18 = $kon_Gpa1GTP_Ste4Ste18 HzPerM\n";
print "koff_Gpa1GTP_Ste4Ste18 = $koff_Gpa1GTP_Ste4Ste18 Hz\n\n";

print "Binding of Ste2-Gpa1-GDP with Ste4-Ste18.\n";
$kon_Ste2Gpa1GDP_Ste4Ste18 = $kon_Gpa1GDP_Ste4Ste18;
$koff_Ste2Gpa1GDP_Ste4Ste18 = $koff_Gpa1GDP_Ste4Ste18 / $Ste2_Gpa1_Ste4Ste18_coop_factor;
print "kon_Ste2Gpa1GDP_Ste4Ste18 = $kon_Ste2Gpa1GDP_Ste4Ste18 HzPerM\n";
print "koff_Ste2Gpa1GDP_Ste4Ste18 = $koff_Ste2Gpa1GDP_Ste4Ste18 Hz\n\n";

print "Binding of Ste2-Gpa1-GTP with Ste4-Ste18.\n";
$kon_Ste2Gpa1GTP_Ste4Ste18 = $kon_Gpa1GDP_Ste4Ste18;
$koff_Ste2Gpa1GTP_Ste4Ste18
    = $GTP_modulating_factor * $koff_Gpa1GDP_Ste4Ste18 
    / $Ste2_Gpa1_Ste4Ste18_coop_factor;
print "kon_Ste2Gpa1GTP_Ste4Ste18 = $kon_Ste2Gpa1GTP_Ste4Ste18 HzPerM\n";
print "koff_Ste2Gpa1GTP_Ste4Ste18 = $koff_Ste2Gpa1GTP_Ste4Ste18 Hz\n\n";
    
print "-----------\n\n";

print "Binding of Ste4 with Ste5.\n\n";

print "Binding of Ste4 with Ste5. (Ste4 not bound to Gpa1.)\n";
$kon_Ste4Ste18_Ste5_HzPerMicromole = 5;
$kon_Ste4Ste18_Ste5 = $kon_Ste4Ste18_Ste5_HzPerMicromole / $molesPerMicromole;
$koff_Ste4Ste18_Ste5 = 2;
print "kon_Ste4Ste18_Ste5 = $kon_Ste4Ste18_Ste5 HzPerM\n";
print "koff_Ste4Ste18_Ste5 = $koff_Ste4Ste18_Ste5 Hz\n\n";

print "Binding of Gpa1-Ste4 with Ste5.\n";
print "My numbers to simulate 'no binding.'\n";
$kon_Gpa1Ste4Ste18_Ste5 = $noBindingOnRate;
$koff_Gpa1Ste4Ste18_Ste5 = $noBindingOffRate;
print "kon_Gpa1Ste4Ste18_Ste5 = $kon_Gpa1Ste4Ste18_Ste5 HzPerM\n";
print "koff_Gpa1Ste4Ste18_Ste5 = $koff_Gpa1Ste4Ste18_Ste5 Hz\n\n";

print "-----------\n\n";

# As supplied
$cell_volume_femtoliters = 50;
$cell_volume = $cell_volume_femtoliters * $litersPerFemtoliter;
print "cell_volume = $cell_volume liters\n\n";

print "-----------\n";
