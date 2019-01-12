# boot_fsc.pl perform parametric bootstrap based on the results of the best run of the best demographic model

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

# Create a bash file to run in parallel fsc in every directory
open (OUT, ">runboot_$ARGV[0].sh") or die;
print OUT "#Run boostrap for $ARGV[0]\n
#PBS -A UQ-SCI-BiolSci\n
#PBS -l select=1:ncpus=9:mem=1GB\n
#PBS -l walltime=1:00:00\n
#PBS -J 1-$ARGV[1]\n
\n
cd /30days/uqharena/bootstrap/$ARGV[0]/$ARGV[0]/$ARGV[0]_\${PBS_ARRAY_INDEX}/; /home/uqharena/fsc26_linux64/fsc26 -t $ARGV[0].tpl -e $ARGV[0].est --initValues $ARGV[0].pv -n100000 -m -L30 -M -c10";

# Launch the job
system ("qsub runboot_$ARGV[0].sh");

close OUT;
exit;

