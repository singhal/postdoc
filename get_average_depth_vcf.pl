use warnings;
use strict;

use List::Util qw(sum);

my $dir = '/ifs/data/c2b2/mp_lab/ss4776/longtailedfinch/masked_genome/';
my $vcf = '/ifs/data/c2b2/mp_lab/ss4776/longtailedfinch/after_vqsr/gatk.ug.ltf.allchrs.snps.indels.vqsr2.vcf';
my $avg_depth = $dir. 'longtailed_avg_depth.txt';
my $summary = $dir . 'longtailed_depth_summary.txt';

sub make_avg_depth {

	#my $call1 = system("gunzip $vcf");
	#$vcf =~ s/\.gz//;

	open(IN, "<$vcf");
	open(OUT, ">$avg_depth");

	my $summed_avg_depth = 0;
	my $num_sites = 0;
	my %cov_hist;	

	while(<IN>) {
		chomp(my $line = $_);
		
		my @d = split(/\t/, $line);

		# don't care about unordered chromosome info
		unless ($d[0] =~ m/Un/ || $d[0] =~ m/random/ || $line =~ m/^#/) {
			# calculate the mean rounded to 3 decimal places for a given base
			
			my $mean = 0;
			if ($line =~ m/DP=([0-9]+)/) {	
				$mean = $1 / 20;
				}

			print OUT $d[0], "\t", $d[1], "\t";
			printf OUT ("%.3f", $mean);
			print OUT "\n"; 
			
			$num_sites++;
			$summed_avg_depth += $mean;
			$cov_hist{int $mean}++;
			}
		}
	close(IN); close(OUT);

	my $cov = $summed_avg_depth / $num_sites;
	open(OUT, ">$summary");
	print OUT "# num sites: $num_sites\n# total_avg_coverage: $summed_avg_depth\n# average_depth: $cov\n";
	print OUT "# histogram of coverage counts\n";
	print OUT "# coverage\tcount\n";
	foreach my $val (sort {$a <=> $b} keys %cov_hist) {
		print OUT $val, "\t", $cov_hist{$val}, "\n";
		} 
	close(OUT);
	}

sub mean {
    return sum(@_)/@_;
}

make_avg_depth();
