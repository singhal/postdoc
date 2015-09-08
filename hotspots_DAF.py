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
parser.add_argument("--sp", help="species for which to run analysis")
args = parser.parse_args()
chr = args.chr
sp = args.sp

genome = '/mnt/gluster/home/sonal.singhal1/reference/ancestral_genome.fa'

hotspots = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/hotspots.LR_seqldhot.ZF_LTF.csv'
out = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/unique_hotspots.%s_%s.mut_skew.csv' % (sp, chr)
if sp == 'LTF':
	vcf = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.%s.coverage.repeatmasked.filtered.phased.vqsr2.vcf.gz' % (chr)
if sp == 'ZF':
	vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.%s.coverage.repeatmasked.filtered.nomendel.shared.noswitch.phased.vqsr2.vcf.gz' % chr

o = open(out, 'w') 
o.write('chr,position,mutation_type,afreq,hotspot_bin\n')

###############

def get_hotspots(hotspots, chr, sp):
	d = pd.read_csv(hotspots)

	if sp == 'ZF':
		lr_high = 'ZF_LR'
		heat_low = 'LTF_heat'
		mid = 'zmid'
	if sp == 'LTF':
		lr_high = 'LTF_LR'
		heat_low = 'ZF_heat'
		mid = 'lmid'

	d = d[d.chr == chr]
	d = d[(d[lr_high] >= 10) & (d[heat_low] <= 1)]

	spots = {}

	for mid in d[mid].dropna().tolist():
		start = mid - 5000
		end = mid + 5000

		for ix, s in enumerate(range(int(start), int(end+1), 100)):
			for bp in range(int(s), int(s+100)):
				spots[bp] = ix
	return spots		


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
spots = get_hotspots(hotspots, chr, sp)

# mutation types
# ancestral is key1, derived is key2
types = {	'A': {'A': 'AT_AT', 'G': 'AT_GC', 'T': 'AT_AT', 'C': 'AT_GC'},
		'C': {'A': 'GC_AT', 'G': 'GC_GC', 'T': 'GC_AT', 'C': 'GC_GC'}, 
		'T': {'A': 'AT_AT', 'G': 'AT_GC', 'T': 'AT_AT', 'C': 'AT_GC'},
		'G': {'A': 'GC_AT', 'G': 'GC_GC', 'T': 'GC_AT', 'C': 'GC_GC'},
	}

v = gzip.open(vcf, 'r')
for l in v:
	if not re.search('#', l):
		d = re.split('\t', l.rstrip())
	
		pos = int(d[1])
		if pos in spots:	
			keep = True
			# check ancestral defined
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
					o.write('%s,%s,%s,%s,%s\n' % (chr, pos, type, derived_af, spots[pos]))
v.close()
o.close()
