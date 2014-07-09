use warnings;
use strict;

my @chrs = qw(1 1A 1B 2 3 4 4A 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 LG2 LG5 LGE22 Z M Un 1A_random 1B_random 1_random 2_random 3_random 4A_random 4_random 5_random 6_random 7_random 8_random 9_random 10_random 11_random 12_random 13_random 14_random 15_random 16_random 17_random 18_random 19_random 20_random 21_random 22_random 23_random 24_random 25_random 26_random 27_random 28_random LGE22_random Z_random);
my $ref_genome = '/KG/sannareddyk/LTFinch/reference/taeGut1.bamorder.fasta';
my $original_bam = '/KG/sannareddyk/LTFinch/realigned_mateFixed_bams/DB.realigned.mateFixed.bam';
my $results_dir = '/mnt/lustre/home/sonal.singhal1/DBF/';
my $gatk = '/mnt/lustre/home/sonal.singhal1/GenomeAnalysisTK.jar';
my $valid_snps = 'XX';

sub generate_initial_variant_files {
	my $out_dir = $results_dir . 'initial_variants/';
        mkdir($out_dir) unless (-d $out_dir);

	foreach my $chr (@chrs) {
		my $out = $out_dir . 'DB.chr' . $chr . '.raw.snps.indels.vcf';
		
		unless(-f $out) {
			my $call = "/home/shyamg/bin/java -Xmx8g -jar $gatk -T UnifiedGenotyper -R $ref_genome -I $original_bam -glm BOTH -L chr" . $chr . " -hets 0.04 -indelHeterozygosity 0.01 -minIndelCnt 1 --output_mode EMIT_VARIANTS_ONLY -o $out -nct 8";
		
		
			my $name = 'initial_variants_chr' . $chr;
			my $out_sh = $out_dir . $name . ".sh";
			open(OUT, ">$out_sh");
			print OUT $call, "\n";
			close(OUT);
			my $call1 = `chmod a+x $out_sh`;
			my $call2 = system("echo $out_sh | qsub -l h_vmem=10g -cwd -V -j y -N \"" . $name . "\"");
			}		
		}
	}

sub concat_variant_files1 {
	my $vcf_dir = $results_dir . 'initial_variants/';
	my $final_file = $vcf_dir . 'DB.allchrs.raw.snps.indels.vcf';
	my $tmp_file = $vcf_dir . 'DB.allchrs.raw.snps.indels.tmp.vcf';
	my $header_file = $vcf_dir . 'header.vcf';

	foreach my $chr (@chrs) {
		my $file =  $vcf_dir . 'DB.chr' . $chr . '.raw.snps.indels.vcf';
		my $call1 = system("grep -v '^#' $file >> $tmp_file");
		my $call2 = system("grep '^#' $file > $header_file") unless(-f $header_file);
		}

	my $final_call = system("cat $header_file $tmp_file > $final_file");
	unlink($tmp_file); unlink($header_file);
	}

sub run_bqsr {
	my $initial_vcf = $results_dir . 'initial_variants/' . 'DB.allchrs.raw.snps.indels.vcf';
	my $out_dir = $results_dir . 'after_bqsr/';
	mkdir($out_dir) unless(-d $out_dir);

	my $recal_grp = $out_dir . 'DB.allchrs.realigned.mateFixed.recal.grp';
	my $out_bam = $out_dir . 'DB.allchrs.realigned.mateFixed.recal.bam';

	my $call1 = "/home/shyamg/bin/java -Xmx8g -jar $gatk -T BaseRecalibrator -R $ref_genome -I $original_bam -knownSites $initial_vcf -o $recal_grp -nct 8";
	my $call2 = "/home/shyamg/bin/java -Xmx8g -jar $gatk -T PrintReads -R $ref_genome -I $original_bam -BQSR $recal_grp -o $out_bam -nct 8";

	my $out_sh = $out_dir . 'run_bqsr.sh';
	open(OUT, ">$out_sh"); print OUT $call1, "\n", $call2, "\n"; close(OUT);
	my $call3 = system("chmod a+x $out_sh");
	my $call4 = system("echo $out_sh | qsub -l h_vmem=10g -cwd -V -j y -N \"run_bqsr\"");
	}

sub generate_bqsr_variant_files {
	my $out_dir = $results_dir . 'after_bqsr/';
	my $recal_bam = $out_dir . 'DB.allchrs.realigned.mateFixed.recal.bam';

        foreach my $chr (@chrs) {                
		my $out = $out_dir . 'DB.chr' . $chr . '.bqsr.snps.indels.vcf';                
		unless(-f $out) {
                        my $call = "/home/shyamg/bin/java -Xmx8g -jar $gatk -T UnifiedGenotyper -R $ref_genome -I $recal_bam -glm BOTH -L chr" . $chr . " -hets 0.04 -indelHeterozygosity 0.01 -minIndelCnt 1 --output_mode EMIT_VARIANTS_ONLY -o $out -nct 9";
                        my $name = 'bqsr_chr' . $chr;                        
			my $out_sh = $out_dir . $name . ".sh";
                        open(OUT, ">$out_sh");
                        print OUT $call, "\n";
                        close(OUT);
                        my $call1 = `chmod a+x $out_sh`;
                        my $call2 = system("echo $out_sh | qsub -l h_vmem=10g -cwd -V -j y -N \"" . $name . "\"");
                        }
                }
        }

sub validate_bqsr_variant_files {
	my $vcf_dir = $results_dir . 'after_bqsr/';

	foreach my $chr (@chrs) {
		my $vcf_file =  $vcf_dir . 'DB.chr' . $chr . '.bqsr.snps.indels.vcf';
		my $name = 'validate_chr' . $chr;
		my $out_sh = $vcf_dir . $name . ".sh";
		
		my $call = "/home/shyamg/bin/java -Xmx4g -jar $gatk -T ValidateVariants -R $ref_genome --variant $vcf_file --validationType NONE";

		open(OUT, ">$out_sh");
		print OUT $call, "\n";
		close(OUT);
		my $call1 = `chmod a+x $out_sh`;
		my $call2 = system("echo $out_sh | qsub -l h_vmem=6g -cwd -V -j y -N \"" . $name . "\"");
		}
	}

sub concat_variant_files2 {
        my $vcf_dir = $results_dir . 'after_bqsr/';
        my $final_file = $vcf_dir . 'DB.allchrs.bqsr.snps.indels.vcf';
        my $tmp_file = $vcf_dir . 'DB.allchrs.bqsr.snps.indels.tmp.vcf';
        my $header_file = $vcf_dir . 'header.vcf';

        foreach my $chr (@chrs) {
                my $file =  $vcf_dir . 'DB.chr' . $chr . '.bqsr.snps.indels.vcf';
                my $call1 = system("grep -v '^#' $file >> $tmp_file");
                my $call2 = system("grep '^#' $file > $header_file") unless(-f $header_file);
                }
        my $final_call = system("cat $header_file $tmp_file > $final_file");
        unlink($tmp_file); unlink($header_file);        
	}

sub run_vqsr {
	my $vcf_dir = $results_dir . 'after_bqsr/';
        my $bqsr_vcf = $vcf_dir . 'DB.allchrs.bqsr.snps.indels.vcf.gz';
	my $out_dir = $results_dir . 'after_vqsr/';
	mkdir($out_dir) unless(-d $out_dir);

	my $trusted_vcf_snp = $results_dir . 'trusted_variants/DB.cortex.raw.snps.vcf.gz';
	my $trusted_vcf_indel = $results_dir . 'trusted_variants/DB.cortex.raw.indels.vcf.gz';
 	my $recalSNP_file = $out_dir . 'DB.allchrs.recal_vqsr.snps.out';
	my $tranchesSNP_file = $out_dir . 'DB.allchrs.tranches_vqsr.snps.out';
	my $rscriptSNP = $out_dir . 'DB.allchrs.vqsr.snps.R';
	my $vqsr_snp = $out_dir . 'DB.allchrs.recal_snps.indels.vcf.gz';
	my $recalIndel_file = $out_dir . 'DB.allchrs.recal_vqsr.indel.out';
	my $tranchesIndel_file = $out_dir . 'DB.allchrs.tranches_vqsr.indel.out';
	my $rscriptIndel = $out_dir . 'DB.allchrs.vqsr.indel.R';
	my $vqsr_snp_indel = $out_dir . 'DB.allchrs.vqsr.snps.indels.vcf';
	
	my $out = $out_dir . 'run_vqsr.sh';
	open(OUT, ">$out");
	print OUT "/home/shyamg/bin/java -Xmx20g -jar $gatk -T VariantRecalibrator -R $ref_genome -input $bqsr_vcf -resource:CORTEX_RAW_INTERSECTION_SNP,known=true,training=true,truth=true,prior=10.0 $trusted_vcf_snp -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an DP -recalFile $recalSNP_file -tranchesFile $tranchesSNP_file -rscriptFile $rscriptSNP -mode SNP -nt 8\n";
	print OUT "/home/shyamg/bin/java -Xmx20g -jar $gatk -T ApplyRecalibration -R $ref_genome -input $bqsr_vcf --ts_filter_level 99.0 -recalFile $recalSNP_file -tranchesFile $tranchesSNP_file -o $vqsr_snp -mode SNP -nt 8\n";
	print OUT "/home/shyamg/bin/java -Xmx20g -jar $gatk -T VariantRecalibrator -R $ref_genome -input $vqsr_snp -resource:CORTEX_RAW_INTERSECTION_INDEL,known=true,training=true,truth=true,prior=10.0 $trusted_vcf_indel -an QD -an MQRankSum -an ReadPosRankSum -an FS -an DP -recalFile $recalIndel_file  -tranchesFile $tranchesIndel_file -rscriptFile $rscriptIndel -mode INDEL -nt 8\n";
	print OUT "/home/shyamg/bin/java -Xmx20g -jar $gatk -T ApplyRecalibration -R $ref_genome -input $vqsr_snp --ts_filter_level 99.0 -recalFile $recalIndel_file -tranchesFile $tranchesIndel_file -o $vqsr_snp_indel -mode INDEL -nt 8\n";
	close(OUT);
	my $call = system("chmod a+x $out");
	}

#generate_initial_variant_files();
#concat_variant_files1();
#run_bqsr();
#generate_bqsr_variant_files();
#validate_bqsr_variant_files();
#concat_variant_files2();
run_vqsr();
