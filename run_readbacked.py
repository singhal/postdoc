import re
import subprocess
import glob 

files = glob.glob('/mnt/lustre/home/sonal.singhal1/DBF/after_vqsr/by_chr/*biallelicSN*gz')

undone = ['chr12', 'chr21', 'chr25', 'chr4', 'chrLG2']

for vcf_file in files:
	chr = re.search('(chr[A-Z|0-9|a-z]+)', vcf_file).group(1)
	if chr in undone:
		out = '/mnt/lustre/home/sonal.singhal1/DBF/phasing/gatk.ug.dbf.%s.phased.vcf.gz' % chr
		call = '/home/shyamg/bin/java -Xmx8g -jar /mnt/lustre/home/sonal.singhal1/bin/GenomeAnalysisTK.jar -T ReadBackedPhasing -R /mnt/lustre/home/sonal.singhal1/reference/taeGut1.bamorder.fasta -I /mnt/lustre/home/sonal.singhal1/DBF/bams/DB.allchrs.realigned.mateFixed.recal.bam --variant %s -L %s -o %s --phaseQualityThresh 20.0' % (vcf_file, vcf_file, out)

		subprocess.call('echo \"%s\" | qsub -l h_vmem=10g -cwd -V -j y -N \'%s\'' % (call, 'rbp_' + chr), shell=True)
