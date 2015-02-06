import re
import os
import subprocess
import glob
import gzip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--sp", help="species for which to run analysis")
parser.add_argument("--chr", help="chromosome for which to run analysis")
args = parser.parse_args()
sp = args.sp
chr = args.chr

ancestral_genome = '/mnt/gluster/home/sonal.singhal1/reference/ancestral_genome.fa'
if chr == 'chrZ':
	if sp == 'ZF':
		vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.chrZ.coverage.repeatmasked.filtered.recodedsex.vqsr.phased.vcf.gz'
	if sp == 'LTF':
		vcf = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.chrZ.filtered.coverage.repeatmasked.recodedsex.vqsr.phased.vcf.gz'
else:
	if sp == 'LTF':
		vcf = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.%s.filtered.coverage.repeatmasked.vqsr.phased.vcf.gz' % chr
	if sp == 'ZF':
		vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.%s.coverage.repeatmasked.filtered.nomendel.phased.vcf.gz' % chr
out = '/mnt/gluster/home/sonal.singhal1/%s/analysis/sfs/%s.sfs' % (sp, chr)

def get_chromosome(genome, chr):
        outfile = genome + '_' + chr
        subprocess.call('~/bin/samtools-0.1.19/samtools faidx %s %s > %s' % (genome, chr, outfile)
, shell=True)
        out_f = open(outfile, 'r')
        chromosome = ''
        locus_name = out_f.next()
        for l in out_f:
                chromosome = chromosome + l.rstrip().upper()
        out_f.close()
        os.remove(outfile)
        return list(chromosome)

o = open(out, 'w')
chr_as_list = get_chromosome(ancestral_genome, chr)
f = gzip.open(vcf, 'r')
for l in f:
	if not re.search('^#', l):
		d = re.split('\t', l.rstrip())
		# only want simple biallelic alleles
		if len(d[3]) == 1 and len(d[4]) == 1:
			genos = []
			for geno in d[9:]:
				geno = re.search('^([^:]+)', geno).group(1)
				genos += re.split('[\||\/]', geno)
			# don't want any missing sites
			if genos.count('.') == 0:
				anc_allele = chr_as_list[int(d[1]) - 1]
				af = 0
				if anc_allele == d[3]:
					af = genos.count('1') / float(len(genos))
				elif anc_allele == d[4]:
					af = genos.count('0') / float(len(genos))
				if af > 0 and af < 1:
					# only want segregating sites
        	                        o.write('%s,%s,%.3f\n' % (chr, d[1], af))
				# don't want any of the other types of sites
f.close()
o.close()
