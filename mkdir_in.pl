# mkdir_in.pl creates directories containing the required input files for launching independent runs per model

# Use: perl mkdir_in.pl (0)number_of_models (1)number_of_independent_runs_per_model

# Assumptions: 
# - The location from where the script is invoked should contain as many directories as models and they should be named as consecutive numbers from 1 to the max number of models. 
# - Each model directory should contain the corresponding SFS (*.obs), estimation (*.est), and template (*.tpl) files.

use warnings;
use strict;

my $dir = 1;
# Create a folder containing the input files for every independent run per model.
while ($dir <= $ARGV[0]) {
	chdir($dir);
	my $r = 1;
	# Create the specified number of folder (independent runs) per model and copy the input files into them.
	while ($r <= $ARGV[1]){
		system "mkdir run$r";
		system "cp *.est *.tpl *.obs run$r";
		$r += 1;
	}
	$dir += 1;
	chdir("..");
}

exit;
