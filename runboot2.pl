# boot_fsc.pl perform parametric bootstrap based on the results of the best run of the best demographic model

# Use: perl boot_fsc.pl (0)root_name (1)number_of_replicates

use warnings;
use strict;

# Extract the parameter values along all the replicates in a new file
open (OUT, ">$ARGV[0]_boot.txt") or die;

my $r = 1;
while ($r <= $ARGV[1]) {
  open (PAR, "$ARGV[0]/$ARGV[0]_$r/$ARGV[0]/$ARGV[0].bestlhoods") or die;
  my $line = <PAR>; # Discard the header line
  $line = <PAR>; # Read the second line that contains the actual parameter values
  chomp $line;
  print OUT "$line\n";
  close PAR;
  $r += 1; # Move to the next replicate
}

close OUT;
exit;
