# extract_ml.pl retrieves the maximum likelihood values of all the independent fastsimcoal runs per model in a table and extracts the key summary files of the best run per model in a new directory

# Use: perl extract_ml.pl (0)number_of_models (1)number_of_independent_runs_per_model (2)root_name_of_files

# Assumptions: 
# - This script should be invoked from the location that contains all the models in separate directories named as consecutive number from 1 to the maximum number of different models after running fastsimcoal.

use warnings;
use strict;

# Create a directory to copy the *.bestlhoods and *.brent_lhoods files of the best run per model
system ("mkdir results_$ARGV[2]");

# Create a file within the above created directory to summarise the likelihood values of all runs 
open (OUT, ">results_$ARGV[2]/$ARGV[2]_bestlhoods.txt") or die;

# Read through all runs and pick up the best run per model while summarising overal results
my $dir = 1;
while ($dir <= $ARGV[0]) { # Do this action for all models
	print OUT "model$dir"; # Print model name in the first column
	my $r = 1;
	my $max = -10000000000000; # Arbitrary very small likelihood value
	my $runML;
	while ($r <= $ARGV[1]){ # Do this action for all runs of a model
		open (L, "$dir/run$r/$ARGV[2]/$ARGV[2].bestlhoods");
		my $line = <L>;
		
    # Extract the maximum likelihood value per run in the output file
		while (my $line = <L>) {
			chomp $line;
			my @col = split(/\t/, $line);
			print OUT "\t$col[-2]"; # The maximum likelihood value is in the column before the last column
			
      # Keep the name and maximum likelihood value of this run if it represents the overall maximum likelihood value so far
      if ($col[-2] > $max) {
				$max = $col[-2];
				$runML = "run$r";
				}
		}
		$r += 1; # Move to the next run
		close L;
	}
	print OUT "\t$runML\n"; # Print the name of the run with the highest maximum likelihood value in the last column
	
  # Copy the *.bestlhoods and *.brent_lhoods files of the best run renaming them to indicate the corresponding model
	system ("cp $dir/$runML/$ARGV[2]/$ARGV[2].bestlhoods results_$ARGV[2]/$ARGV[2]_model$dir.bestlhoods");
	system ("cp $dir/$runML/$ARGV[2]/$ARGV[2].brent_lhoods results_$ARGV[2]/$ARGV[2]_model$dir.brent_lhoods");
	
  $dir += 1; # Move to the next model
}

close OUT;
exit;
