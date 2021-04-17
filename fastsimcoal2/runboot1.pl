# runboot1.pl simulates SFS files based on an given demographic model for performing parametric bootstrap

# Use: perl boot_fsc.pl (0)root_name (1)number_of_replicates

# Assumptions: 
# - The parameter (*.par), initial parameter values (*.pv), template (*.tpl), and estimation (*.est) files should be contained in the local directory

use warnings;
use strict;

# Run fsc to simulate the specified number of SFS from a parameter file
system ("[PATH]/fsc26 -i $ARGV[0].par -n$ARGV[1] -j -m -s0 -x -I -q");

# Copy the initial parameter values, template, and estimation files to the newly created directories containing the simulated SFS
my $rep = 1;
while ($rep <= $ARGV[1]) {
  system ("cp $ARGV[0].tpl $ARGV[0]/$ARGV[0]_$rep/");
  system ("cp $ARGV[0].est $ARGV[0]/$ARGV[0]_$rep/");
  system ("cp $ARGV[0].pv $ARGV[0]/$ARGV[0]_$rep/");
  $rep += 1; # Move to the next replicate directory
}

exit;

