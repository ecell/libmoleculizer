#!/usr/bin/perl

($dmp_header) = @ARGV;

my $total = 0.0;
my $totSq = 0.0;
my $count = 0;

while(my $value = <STDIN>)
{
  chomp $value;

  $total += $value;
  $totSq += ($value * $value);
  $count++;
}

my $mean = $total / $count;
my $second_moment = $totSq / $count;
my $sampleVariance
  = $count * ($second_moment - ($mean * $mean)) / ($count - 1);

print "$dmp_header\t$mean\t$sampleVariance\n";
