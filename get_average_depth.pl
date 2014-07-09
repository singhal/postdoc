use warnings;
use strict;

use List::Util qw(sum);

my $dir = '/mnt/lustre/home/sonal.singhal1/DBF/masked_genome/';
my $depth = $dir . 'doublebarred_depth.txt';
my $bam_list = $dir . 'doublebarred_bam_list.txt';
my $avg_depth = $dir. 'doublebarred_avg_depth.txt';
my $summary = $dir . 'doublebarred_depth_summary.txt';

sub make_sample_depth {
	}

sub make_avg_depth {
	open(IN, "<$depth");
	open(OUT, ">$avg_depth");

	my $summed_avg_depth = 0;
	my $num_sites = 0;
	my %cov_hist;	

	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/, $line);
		
		# don't care about unordered chromosome info
		unless ($d[0] =~ m/Un/ || $d[0] =~ m/random/) {
			# calculate the mean rounded to 3 decimal places for a given base
			# logic ignores the first two fields because they indicate the chromosome and base pos
			my $mean = 0;
			if (scalar(@d[2 .. $#d]) > 0) {
				$mean = mean(@d[2 .. $#d]);
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
