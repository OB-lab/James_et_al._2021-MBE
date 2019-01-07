# getpopmap.pl generates a tab-delimited population specification file from a CVF file containing either two or three populations

# Use: perl getpopmap.pl (0)input_vcf_file (1)output_popmap_file (2)name_pop1 (3)name_pop2 (4)name_pop3

# Assumptions: 
# - Sample names are composed by PopName-SampleID_PopName-SampleID (i.e. D12-291_D12-291). Most importantly, 
#   the first part of the sample name should contain the population name followed by a hyphen.
# - The header line containing the sample names is tab-delimitated and stats by #CHROM
# - Sample names start at column number 10

use warnings;
use strict;

open (VCF, $ARGV[0]) or die; # Load the input VCF file
open (OUT, ">$ARGV[1]") or die; # Create the output file

# Read the VCF file line per line until reaching the header and extracting the sample and population names
my $c=0; # Keep reading the file until finding the header
while ($c==0) {
	my $l = <VCF>; # Read the line
	chomp $l;
	my @col = split(/\t/, $l);
	if ($col[0] eq "#CHROM") { # Check whether the current line corresponds to the header
		for (my $i=9; $i<scalar @col; $i++) { # Read sample names from header's column number 10 onwards
			my @sn = split (/-/, $col[$i]); # Extract population name from the sample name
			# If there are two populations, write the sample and population names in a new line in the output file...
			if (scalar(@ARGV) == 4) {
				if ($sn[0] eq $ARGV[2]) { # if the sample belongs to population 1
					print OUT "$col[$i]\t$ARGV[2]\n"; 
				}
				if ($sn[0] eq $ARGV[3]) { # if the sample belongs to population 2
					print OUT "$col[$i]\t$ARGV[3]\n";
				}
			}
			# If there are three populations, write the sample and population names in a new line in the output file...
			if (scalar(@ARGV) == 5) { 
				if ($sn[0] eq $ARGV[2]) { # if the sample belongs to population 1
					print OUT "$col[$i]\t$ARGV[2]\n"; 
				}
				if ($sn[0] eq $ARGV[3]) { # if the sample belongs to population 2
					print OUT "$col[$i]\t$ARGV[3]\n";
				}
				if ($sn[0] eq $ARGV[4]) { # if the sample belongs to population 3
					print OUT "$col[$i]\t$ARGV[4]\n";
				}
			}	
		}
		$c = 1; # Stop reading the file because the header's useful information was already taken
		}
}

close VCF;
close OUT;
exit;
