# getpopmap.pl
# 27/11/2018
# henry.arenasc@gmail.com
# Generte tab-delimited population map file from a CVF file.
# Use: per getpopmap.pl (0)vcf_file (1)output_filename (2)root_name_pop1 (3)root_name_pop2

use warnings;
use strict;

open (CVF, $ARGV[0]) or die;
open (OUT, ">$ARGV[1]") or die;

my $c=0;
while ($c==0) {
	my $l = <CVF>;
	chomp $l;
	my @col = split(/\t/, $l);
	if ($col[0] eq "#CHROM") {
		for (my $i=9; $i<scalar @col; $i++) { #Assuming that pop names start at column #10 in the VCF file
			my @sn = split (/-/, $col[$i]);
			if ($sn[0] eq $ARGV[2]) {
				print OUT "$col[$i]\t$ARGV[2]\n";
			}
			if ($sn[0] eq $ARGV[3]) {
				print OUT "$col[$i]\t$ARGV[3]\n";
			}
		}
		$c = 1;
		}
}

close CVF;
close OUT;
exit;
