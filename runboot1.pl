# runboot1.pl perform parametric bootstrap based on the results of the best run of the best demographic model

# Use: perl boot_fsc.pl (0)root_name (1)number_of_replicates

use warnings;
use strict;

# Run fsc to simulate the specified number of SFS from the *.par file
system ("/home/uqharena/fsc26_linux64/fsc26 -i $ARGV[0].par -n$ARGV[1] -j -m -s0 -x -I -q");

# Copy the *.tpl, *.est, and *.pv files to the replicate directories containing the simulated SFS
my $rep = 1;
while ($rep <= $ARGV[1]) {
  system ("cp $ARGV[0].tpl $ARGV[0]/$ARGV[0]_$rep/");
  system ("cp $ARGV[0].est $ARGV[0]/$ARGV[0]_$rep/");
  system ("cp $ARGV[0].pv $ARGV[0]/$ARGV[0]_$rep/");
  $rep += 1; # Move to the next replicate directory
}

exit;

