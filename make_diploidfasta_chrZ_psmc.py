from itertools import izip
import gzip
import re
import sys
import subprocess
import os
import argparse


def get_vcf(vcffile):
	code = {'A': {'A': 'A', 'C': 'M', 'G': 'R', 'T': 'W'},
		'C': {'A': 'M', 'C': 'C', 'G': 'S', 'T': 'Y'},
		'G': {'A': 'R', 'C': 'S', 'G': 'G', 'T': 'K'},
		'T': {'A': 'W', 'C': 'Y', 'G': 'K', 'T': 'T'}}

        file = gzip.open(vcffile, 'r')
        var = dict()
	ids = []
        for l in file:
                if re.match('#CHROM', l):
			d = re.split('\t', l.rstrip())
			ids = d[9:len(d)]
			for id in ids:
				var[id] = dict()
		if not re.match("#", l):
                        d = re.split('\t', l.rstrip())
			all_num = dict()
			alleles = [d[3], d[4]]
			for ix, allele in enumerate(alleles):
				all_num[str(ix)] = allele
			# populate var
			for id, geno in zip(ids, d[9:len(d)]):
				geno1 = re.search('(\S)\/', geno).group(1)
				geno2 = re.search('\/(\S)', geno).group(1)
				if geno1 != '.':
					var[id][int(d[1]) - 1] = code[all_num[geno1]][all_num[geno2]]
				else:
					var[id][int(d[1]) - 1] = 'N'
        file.close()
        return var


def get_chromosome(genome, chr):
        outfile = genome + '_' + chr
        subprocess.call('~/bin/samtools-0.1.19/samtools faidx %s %s > %s' % (genome, chr, outfile), shell=True)
        out_f = open(outfile, 'r')
        chromosome = ''
        locus_name = out_f.next()
        for l in out_f:
                chromosome = chromosome + l.rstrip()
        out_f.close()
        os.remove(outfile)
        return list(chromosome)


def print_seq(var, chr_as_list, masked, out_dir, chr):
        for ind in var:
		out_file = '%s%s_chrZ.fasta' % (out_dir, ind)
		out_f = open(out_file, 'a')
                out_f.write('>%s\n' % chr)
                tmp_chr = list(chr_as_list)
                for pos, base in enumerate(masked):
                        if int(base) > 0:
                                tmp_chr[pos] = 'N'
                for pos in var[ind]:
                        tmp_chr[pos] = var[ind][pos]
                for i in xrange(0, len(tmp_chr), 60):
                        out_f.write(''.join(tmp_chr[i:i+60]) + '\n')
        out_f.close()
	return


def main():
	sp = 'LTF'
	out_dir = '/mnt/gluster/home/sonal.singhal1/%s/analysis/PSMC/' % sp
	if sp == 'ZF':
		vcf_in = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/all_vcf/for_shapeit/gatk.ug.finch19.chrZ.allfilters.recodedsex.recoded_biallelicSNPs.males.vcf.gz'
		masked_genome = '/mnt/gluster/home/sonal.singhal1/ZF/masked_genome/ZF.masked_genome.repeat_masked.switch_masked.fa'
	elif sp == 'LTF':
		vcf_in = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/for_shapeit/gatk.ug.ltf.chrZ.allfilters.recodedsex.recoded_biallelicSNPs.males.vcf.gz'
		masked_genome = '/mnt/gluster/home/sonal.singhal1/LTF/masked_genome/LTF.masked_genome.repeat_masked.fa'
	genome_ref = '/mnt/gluster/home/sonal.singhal1/reference/taeGut1_60.bamorder.fasta'
	
	var = get_vcf(vcf_in)
	masked = get_chromosome(masked_genome, 'chrZ')
	chr_ref = get_chromosome(genome_ref, 'chrZ')
	print_seq(var, chr_ref, masked, out_dir, 'chrZ')
		
if __name__ == "__main__":
    main()
