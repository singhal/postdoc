import pandas as pd
import re
import numpy as np
from itertools import izip
import argparse
import gzip
import os
import sys
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("--chr", help="chromosome for which to run analysis")
args = parser.parse_args()
chr = args.chr

sp = 'LTF'

genome = '/mnt/gluster/home/sonal.singhal1/reference/ancestral_genome.fa'
gff = '/mnt/gluster/home/sonal.singhal1/reference/Taeniopygia_guttata.gff'
# long autosomal chromosomes
file = '/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/maps/%s_recombination_bpen100.txt' % (sp, chr)
out = '/mnt/gluster/home/sonal.singhal1/%s/analysis/mut_skew/%s.recombination_mut_skew.csv' % (sp, chr)
if sp == 'LTF':
	vcf = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.%s.coverage.repeatmasked.filtered.phased.vqsr2.vcf.gz' % (chr)
if sp == 'ZF':
	vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.%s.coverage.repeatmasked.filtered.nomendel.shared.noswitch.phased.vqsr2.vcf.gz' % chr

o = open(out, 'w') 
o.write('chr,position,mutation_type,afreq,distance_to_TSS,rho\n')

###############

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

# ancestral genome
chr_as_list = get_chromosome(genome, chr)

# mutation types
# ancestral is key1, derived is key2
types = {	'A': {'A': 'AT_AT', 'G': 'AT_GC', 'T': 'AT_AT', 'C': 'AT_GC'},
		'C': {'A': 'GC_AT', 'G': 'GC_GC', 'T': 'GC_AT', 'C': 'GC_GC'}, 
		'T': {'A': 'AT_AT', 'G': 'AT_GC', 'T': 'AT_AT', 'C': 'AT_GC'},
		'G': {'A': 'GC_AT', 'G': 'GC_GC', 'T': 'GC_AT', 'C': 'GC_GC'},
	}

# get in genes
d = pd.read_csv(gff, sep='\t', header=0, names=['chr', 'type', 'cds_mrna', 
                                                'start', 'stop', 'score', 
                                                'orientation', 'codon_pos', 'id'])
chr_short = chr.replace('chr','')
d = d[d.chr == chr_short]
gdata = d[ d.cds_mrna == 'mRNA' ].set_index('id').to_dict()
genes = d[ d.cds_mrna == 'mRNA' ].id.unique()
tss = {}
for gene in genes:
        orientation = 1
        if gdata['orientation'][gene] == '+':
                tss_start = gdata['start'][gene]
        elif gdata['orientation'][gene]  == '-':
                tss_start = gdata['stop'][gene]
                orientation = -1
        tss[ tss_start ]  = orientation
tss_starts = sorted(tss.keys())

# rho
rhos = {}
rho_df = pd.read_csv(file, sep=" ", skiprows=3, header=None,
        names=['left_snp', 'right_snp', 'meanrho', 'p0.025', 'p0.975'])
for start, rho in izip(rho_df.left_snp, rho_df.meanrho):
        rhos[start] = rho
rho_starts = [None] + sorted(rhos.keys())

tss_ix = 0
rho_ix = 0

v = gzip.open(vcf, 'r')
for l in v:
	if not re.search('#', l):
		d = re.split('\t', l.rstrip())
		
		keep = True
		
		# check ancestral defined
		pos = int(d[1])
		ancestral = 'N'
		if chr_as_list[pos -1] != 'N':
			ancestral = chr_as_list[pos -1] 
		else:
			keep = False

		# check not indel
		alleles = [d[3]] + re.split(',', d[4])
		for allele in alleles:
			if len(allele) > 1:
				keep = False

		# check number of alleles
		genos = []
		for geno in d[9:]:
			geno = re.search('^([^:]+)', geno).group(1)
			genos += re.split('[|/]', geno)
		genos  = [x for x in genos if not re.search('\.', x)]

		af = {}
		for ix, allele in enumerate(alleles):
			if genos.count(str(ix)) > 0:
				af[ allele ] = genos.count(str(ix)) / float(len(genos))
		# only want monoallelics or biallelic
		if len(af) > 2:
			keep = False

		# now generate the hash
		if keep:
			derived_af = None
			type = None
			if ancestral in af:
				derived_af = 1 - af[ancestral]
				derived = [x for x in af.keys() if x != ancestral]
				if len(derived) > 0:
					type = types[ancestral][derived[0]]
			else:
				if len(af) == 1:
					derived_af = 1
					type = types[ancestral][af.keys()[0]]

			if type:
				if (tss_ix + 1) < len(tss_starts):
					midpoint = int((tss_starts[tss_ix] +  tss_starts[tss_ix + 1]) / 2.0)
					if pos >= midpoint:
						tss_ix += 1
				if (rho_ix + 1) < len(rho_starts):
					if pos >= rho_starts[rho_ix + 1]:
						rho_ix += 1

				# get the min distance to tss with orientation
				tss_start = tss_starts[tss_ix]
				min_dist = (pos - tss_start) * tss[ tss_start ]

				# get the rho
				if rho_ix > 0:
					pos_rho = rhos[rho_starts[rho_ix]]
					o.write('%s,%s,%s,%s,%s,%s\n' % (chr, pos, type, derived_af, min_dist, pos_rho))
v.close()
o.close()


