import re
import glob
import sys
import subprocess
import os
import copy
import argparse
import gzip
import random

def parse_haps(hap_file, chr):
	f = gzip.open(hap_file, 'r')
	
	var = {0: dict(), 1: dict()}

	for l in f:
		if not re.search('^#', l):
			d = re.split('\t', l.rstrip())
			alleles = [d[3]] + re.split(',', d[4])
			indel = False
			for allele in alleles:
				if len(allele) > 1:
					indel = True
			if indel is False:
				alldict = {}
				for ix, allele in enumerate(alleles):
					alldict[ix] = allele
	
				geno1 = alldict[ int(re.search('(\d)\/', l).group(1)) ]
				geno2 = alldict[ int(re.search('\/(\d)', l).group(1)) ]

				pos = int(d[1]) - 1

				if random.uniform(0,1) > 0.5:
					var[0][pos] = geno1
					var[1][pos] = geno2
				else:
					var[0][pos] = geno2
					var[1][pos] = geno1
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
		for pos, base in enumerate(masked):
			if base in ['4', '5', '6', '7']:
				tmp_chr[pos] = 'N'
		for pos in var[ind]:
                        tmp_chr[pos] = var[ind][pos]
		for i in xrange(0, len(tmp_chr), 60):
        		out_f.write(''.join(tmp_chr[i:i+60]) + '\n')
	out_f.close()


def main():
        parser = argparse.ArgumentParser()
        parser.add_argument("--chr", help="chromosome for which to run analysis")
        args = parser.parse_args()
        chr = args.chr

	hap_file = '/mnt/gluster/home/sonal.singhal1/DBF/after_vqsr/by_chr/gatk.ug.dbf.%s.filtered.coverage.vqsr.vcf.gz' % chr
	out_file = '/mnt/gluster/home/sonal.singhal1/gene_trees/DBF_%s_haplotypes.fasta' % chr
	genome = '/mnt/gluster/home/sonal.singhal1/reference/taeGut1.bamorder.fasta'
	masked_genome = '/mnt/gluster/home/sonal.singhal1/DBF/masked_genome/DBF.masked_genome.fa'
	
	chr_as_list = get_chromosome(genome, chr)
	masked = get_chromosome(masked_genome, chr)
	var = parse_haps(hap_file, chr)
	print_seq(var, chr_as_list, masked, out_file)

if __name__ == "__main__":
    main()
		
