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

			if d[4] != '.' and d[6] == 'PASS':
	                        alleles = [d[3]] + re.split(",", d[4])
                        
                        	# do not want to include indels in this variant dictionary
                        	# cannot use these to polarize
                        	indel = False
                        	for allele in alleles:
                                	if len(allele) > 1:
                                        	indel = True
                        
                        	if not indel:
					all_num = dict()
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
		out_file = '%s%s.fasta' % (out_dir, ind)
		out_f = open(out_file, 'a')
                out_f.write('>%s\n' % chr)
                tmp_chr = list(chr_as_list)
                for pos, base in enumerate(masked):
                        if base in ['4', '5', '6', '7']:
                                tmp_chr[pos] = 'N'
                for pos in var[ind]:
                        tmp_chr[pos] = var[ind][pos]
                for i in xrange(0, len(tmp_chr), 60):
                        out_f.write(''.join(tmp_chr[i:i+60]) + '\n')
        out_f.close()
	return


def main():
	chrs = [ 'chr1', 'chr1A', 'chr1B', 'chr2', 'chr3',  'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
	         'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
	         'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28',
	         'chrLG2', 'chrLG5', 'chrLGE22' ]

	for chr in chrs:
		out_dir = '/mnt/gluster/home/sonal.singhal1/DBF/analysis/PSMC/'
		vcf_in = '/mnt/gluster/home/sonal.singhal1/DBF/after_vqsr/by_chr/gatk.ug.dbf.%s.filtered.coverage.biallelicSNPs.vqsr.vcf.gz' % chr
		masked_genome = '/mnt/gluster/home/sonal.singhal1/DBF/masked_genome/DBF.masked_genome.fa'
		genome_ref = '/mnt/gluster/home/sonal.singhal1/reference/taeGut1_60.bamorder.fasta'
	
		var = get_vcf(vcf_in)
		masked = get_chromosome(masked_genome, chr)
		chr_ref = get_chromosome(genome_ref, chr)
		print_seq(var, chr_ref, masked, out_dir, chr)
		del var
		
if __name__ == "__main__":
    main()
