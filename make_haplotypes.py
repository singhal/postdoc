import re
import glob
import sys
import subprocess
import os
import copy
import argparse
import gzip
import time
import random

def parse_haps(hap_file, chr):
	# check if indel
	# turn numbers into letters
	# check if phased or not
	f = gzip.open(hap_file, 'r')
	
	var = {}
	
	for l in f:
		if not re.search('#', l):
			d = re.split('\t', l.rstrip())
			
			indel = False
			# define a dict with allele number pointing to allele identity
			alleles = [d[3]] + re.split(',', d[4])
			num_allele = {}
			for ix, allele in enumerate(alleles):
				num_allele[str(ix)] = allele
				if len(allele) > 1:
					indel = True
			# for missing data
			num_allele['.'] = 'N'		

			# only working with snps
			if not indel:
				# to keep track of individuals
				tracker = 0
				for geno in d[9:]:
					geno = re.search('^([^:]+)', geno).group(1)
					genos = []
					if re.search('\|', geno):
						genos = re.split('\|', geno)
					elif re.search('\/', geno):
						# for unknown phase, randomize!
						genos = re.split('\/', geno)
						random.shuffle(genos)
					else:
						# support for xchromosomes
						genos = [geno]
					if tracker not in var:
						var[tracker] = {}
					if len(genos) == 2:
						var[tracker][int(d[1])-1] = num_allele[genos[0]]
						if (tracker+1) not in var:
							var[tracker+ 1] = {}
						var[tracker + 1][int(d[1])-1] = num_allele[genos[1]]
					else:
						var[tracker][int(d[1])-1] = num_allele[genos[0]]
					# tracks the haplotype count based on the ploidy of individuals
					tracker += len(genos)
	f.close()

	return var

				
def get_chromosome(genome, chr):
        outfile = genome + '_' + chr
        subprocess.call('~/bin/samtools-0.1.19/samtools faidx %s %s > %s' % (genome, chr, outfile), shell=True)
        out_f = open(outfile, 'r')
	chromosome = ''
        locus_name = out_f.next()
	for l in out_f:
		chromosome = chromosome + l.rstrip().upper()
	out_f.close()
	os.remove(outfile)
        return list(chromosome)


def print_seq(var, chr_as_list, masked, out_file):
	out_f = open(out_file, 'w')
	for ind in var:
		out_f.write('>haplo%s\n' % ind)
		tmp_chr = list(chr_as_list)
		for pos in var[ind]:
			tmp_chr[pos] = var[ind][pos]
		for pos, base in enumerate(masked):
			if base in ['4', '5', '6', '7', '8']:
				tmp_chr[pos] = 'N'
		for i in xrange(0, len(tmp_chr), 60):
        		out_f.write(''.join(tmp_chr[i:i+60]) + '\n')
	out_f.close()


def main():
        parser = argparse.ArgumentParser()
        parser.add_argument("--chr", help="chromosome for which to run analysis")
	parser.add_argument("--sp", help="chromosome for which to run analysis")
        args = parser.parse_args()
        chr = args.chr
	sp = args.sp

	if sp == 'ZF':
		hap_file = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.%s.coverage.repeatmasked.filtered.nomendel.phased.vcf.gz' % chr
		if chr == 'chrZ':
			hap_file = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.chrZ.coverage.repeatmasked.filtered.recodedsex.vqsr.phased.vcf.gz'
	if sp == 'LTF':
		hap_file = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.%s.filtered.coverage.repeatmasked.vqsr.phased.vcf.gz' % chr
		if chr == 'chrZ':
			hap_file = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.chrZ.filtered.coverage.repeatmasked.recodedsex.vqsr.phased.vcf.gz'
	out_file = '/mnt/gluster/home/sonal.singhal1/gene_trees/%s_%s_haplotypes.fasta' % (sp, chr)
	genome = '/mnt/gluster/home/sonal.singhal1/reference/taeGut1.bamorder.fasta'
	masked_genome = '/mnt/gluster/home/sonal.singhal1/%s/masked_genome/%s.masked_genome.repeat_masked.fa' % (sp, sp)
	
	chr_as_list = get_chromosome(genome, chr)
	masked = get_chromosome(masked_genome, chr)
	var = parse_haps(hap_file, chr)
	print_seq(var, chr_as_list, masked, out_file)

if __name__ == "__main__":
    main()
		
