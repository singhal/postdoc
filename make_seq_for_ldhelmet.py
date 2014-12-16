import re
import glob
import sys
import subprocess
import os
import copy
import argparse
import gzip
import time

def get_vcf_var(vcf_file, min_allele):
        site_freq = {}
        vcf_f = gzip.open(vcf_file, 'r')
        
        for l in vcf_f:
                if not re.search('^#', l):
                        d = re.split('\t', l)

                        allele1 = len(re.findall('0\/', l)) + len(re.findall('\/0', l))
                        allele2 = len(re.findall('1\/', l)) + len(re.findall('\/1', l))
                        min_ac = min(allele1, allele2)
			if min_ac <= min_allele:
                        	site_freq[int(d[1])] = 1
        vcf_f.close()

        return site_freq


def parse_haps(hap_file, site_freq, sites_file, chr):
	f = open(hap_file, 'r')
	
	var = {}
	sites = {}

	for l in f:
		d = re.split('\s+', l.rstrip())
		if d[3] in ['A', 'T', 'C', 'G'] and d[4] in ['A', 'T', 'C', 'G']:
			if int(d[2]) not in site_freq:
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
	
	sites_f = open(sites_file, 'w')
	for site in sorted(var[0].keys()):
		sites_f.write('%s,%s\n' % (chr, site + 1))
	sites_f.close()

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
			if int(base) > 3:
				tmp_chr[pos] = 'N'
		for i in xrange(0, len(tmp_chr), 60):
        		out_f.write(''.join(tmp_chr[i:i+60]) + '\n')
	out_f.close()


def main():
        parser = argparse.ArgumentParser()
        parser.add_argument("--chr", help="chromosome for which to run analysis")
        args = parser.parse_args()
        chr = args.chr

	vcf_file = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/for_shapeit/gatk.ug.unrel_zf.%s.coverage.repeatmasked.filtered.recoded_biallelicSNPs.nomendel.vcf.gz' % chr
	hap_file = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/results/%s_haplotypes.haps' % chr
	out_file = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/%s_haplotypes.fasta' % chr
	site_file = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/%s_sites.csv' % chr
	genome = '/mnt/gluster/home/sonal.singhal1/reference/taeGut1.bamorder.fasta'
	masked_genome = '/mnt/gluster/home/sonal.singhal1/ZF/masked_genome/ZF.masked_genome.repeat_masked.fa'
	min_allele = 1

	site_freq = get_vcf_var(vcf_file, min_allele)
	chr_as_list = get_chromosome(genome, chr)
	masked = get_chromosome(masked_genome, chr)
	var = parse_haps(hap_file, site_freq, site_file, chr)
	print_seq(var, chr_as_list, masked, out_file)

if __name__ == "__main__":
    main()
		
