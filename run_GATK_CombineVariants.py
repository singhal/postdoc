import glob
import re
import subprocess

chrs = [ 'chr1', 'chr1A', 'chr1B', 'chr2', 'chr3', 'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
         'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
         'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28',
         'chrLG2', 'chrLG5', 'chrLGE22', 'chrZ' ]

for chr in chrs:
	vcf1 = '/mnt/lustre/home/sonal.singhal1/ZF/phasing/family_approach/family_vcf/by_chr/gatk.ug.zf_fam.%s.coverage.no_mendel.poly_in_any_zf.vqsr.vcf.gz' % chr
	vcf2 = '/mnt/lustre/home/sonal.singhal1/ZF/phasing/family_approach/unrel_vcf/by_chr/gatk.ug.zf_unrel.%s.coverage.no_mendel.poly_in_any_zf.vqsr.vcf.gz' % chr
	out = '/mnt/lustre/home/sonal.singhal1/ZF/phasing/family_approach/gatk.ug.all_zf.%s.coverage.no_mendel.vqsr.vcf.gz' % chr 

	call = '/home/shyamg/bin/java -Xmx5g -jar ~/bin/GenomeAnalysisTK.jar -R ~/reference/Zfinch.fa -T CombineVariants --variant %s --variant %s -o %s -genotypeMergeOptions UNSORTED' % (vcf1, vcf2, out) 

	if chr == 'chr22':
		subprocess.call('echo \"%s\" | qsub -l h_vmem=7g -cwd -V -j y -N \'%s\'' % (call, chr + 'cat_gatk'), shell=True)
