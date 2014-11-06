import re
import glob
import sys
import subprocess
import os
import copy
import argparse
import gzip
import time

def parse_haps(hap_file, chr):
	f = open(hap_file, 'r')
	
	var = {}
	sites = {}

	for l in f:
		d = re.split('\s+', l.rstrip())
		if d[3] in ['A', 'T', 'C', 'G'] and d[4] in ['A', 'T', 'C', 'G']:
			if int(d[2]) in sites:
				var.pop(int(d[2]), None)
			else:
				sites[int(d[2])] = 1
				for ix, i in enumerate(d[5:len(d)]):
					if ix not in var:
						var[ix] = dict()
					if i == '0':
						var[ix][int(d[2]) - 1] = d[3]
					else:
						var[ix][int(d[2]) - 1] = d[4]
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
			if base in ['4', '5', '6', '7', '8']:
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

	hap_file = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/results/%s_haplotypes.haps' % chr
	out_file = '/mnt/gluster/home/sonal.singhal1/gene_trees/ZF_%s_haplotypes.fasta' % chr
	genome = '/mnt/gluster/home/sonal.singhal1/reference/taeGut1.bamorder.fasta'
	masked_genome = '/mnt/gluster/home/sonal.singhal1/ZF/masked_genome/ZF.masked_genome.repeat_masked.fa'
	
	chr_as_list = get_chromosome(genome, chr)
	masked = get_chromosome(masked_genome, chr)
	var = parse_haps(hap_file, chr)
	print_seq(var, chr_as_list, masked, out_file)

if __name__ == "__main__":
    main()
		
