use warnings;
use strict;

my $dir = '/ifs/data/c2b2/mp_lab/ss4776/zebrafinch/vqsr_profiling/';
# my $vqsr_file = '/mnt/lustre/home/sonal.singhal1/DBF/after_vqsr/DB.allchrs.vqsr.snps.indels.vcf.gz';
my $var_file = $dir . 'gatk.ug.unrelzf.allchrs.var.vqsr2.vcf';
my $out = $dir . 'ZF.allchrs.var.vqsr_score_details.csv';

sub make_filtered {
	# my $call = system("zcat $vqsr_file | grep 'AF=' > $var_file");
	}

sub process_var_file {
	open(IN, "<$var_file");
	open(OUT, ">$out");

	my @vals = qw(AC AF AN BaseQRankSum DP Dels FS HaplotypeScore InbreedingCoeff MLEAC MLEAF MQ MQ0 MQRankSum QD ReadPosRankSum VQSLOD);

	print OUT "type,filter,", join(",",@vals), "\n";
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/, $line);
		my $type = 'SNP';
		if ($d[3] =~ m/[ATGC][ATGC]/ || $d[4] =~ m/[ATGC][ATGC]/) {
			$type = 'INDEL';
			}
		my $pass = 'FAIL';
		$pass = 'PASS' if $line =~ m/PASS/;

		my %vals;
		while($d[7] =~ m/([^;]+)=([-|\d|\.]+)/g) {
			$vals{$1} = $2;
			}

		print OUT $type, ",", $pass, ",";
		foreach my $val (@vals) {
			if (exists($vals{$val})){
				print OUT $vals{$val}, ",";
				}
			else {
				print OUT ",";
				}
			}
		print OUT "\n";
		}

	close(IN); close(OUT);
	}

#make_filtered();
process_var_file();
