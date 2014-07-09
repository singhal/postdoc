use warnings;
use strict;

my $vcf1 = '/mnt/lustre/home/sonal.singhal1/DBF/trusted_variants/DB.cortex.snps.indels.vcf';
my $vcf2 = '/mnt/lustre/home/sonal.singhal1/DBF/initial_variants/DB.allchrs.raw.snps.indels.vcf';

my $out = "intersected.vcf";

my %snps;
open(IN1, "<$vcf1");
while(<IN1>) {
	chomp(my $line = $_);
	unless ($line =~ m/^#/) {
		my @d = split(/\t/, $line);
		$snps{$d[0]}{$d[1]} = $d[4];
		}
	}
close(IN1);

print "vcf1 done\n";

open(OUT, ">$out");	
open(IN2, "<$vcf2");
while(<IN2>) {
	chomp(my $line = $_);
	if ($line =~ m/^#/) {
		print OUT $line, "\n";
		}
	else {
		my @d = split(/\t/, $line);
		if (exists($snps{$d[0]}{$d[1]})) {
			if ($snps{$d[0]}{$d[1]} =~ m/\,/ || $d[4] =~ m/\,/) {
				my @snps1 = split(/,/, $snps{$d[0]}{$d[1]});
				my @snps2 = split(/,/, $d[4]);

				my $match = 0;
				foreach my $snp1 (@snps1) {
					foreach my $snp2 (@snps2) {
						$match = 1 if $snp1 eq $snp2;
						}
					}
				print OUT $line, "\n" if $match;
				}
			else {
				if ($snps{$d[0]}{$d[1]} eq $d[4]) {
					print OUT $line, "\n";
					}
				}

			}
		}	
	}
close(IN2);
close(OUT);
