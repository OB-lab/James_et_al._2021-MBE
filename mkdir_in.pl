# mkdir_in.pl
# 12/09/2018
# henry.arenasc@gmail.com
# Create N folder files for equal number of independent runs for a given X number of models
# Use: perl mkdir_in.pl (0) number of models (1) number of independent runs per model

# The *.obs, *.est, and *.tpl files should be in every model's folder named from from 1 to $ARGV[0]

use warnings;
use strict;

my $dir = 1;
# Create a folder containing the input files for every independent run per model.
# Model folders should be named from 1 to n.
while ($dir <= $ARGV[0]) {
	chdir($dir);
	my $r = 1;
	# Create the specified number of folder (independent runs) per model and move the input data into them.
	while ($r <= $ARGV[1]){
		system "mkdir run$r";
		system "cp *.est *.tpl *.obs run$r";
		$r += 1;
	}
	$dir += 1;
	chdir("..");
}

exit;
