import re
import subprocess
from itertools import izip
import pandas as pd
import numpy as np

def get_uniq_spots(d, sp):
	opp = {'ZF': 'LTF', 'LTF': 'ZF'}

	d = pd.read_csv(d)
	
	d = d.dropna()
	d[['ZF_LR', 'LTF_LR']] = d[['ZF_LR', 'LTF_LR']].astype(float)

	d = d[d['%s_LR' % sp] >= 10]
	d = d[d['%s_heat' % opp[sp]] < 5]
	d = d[d['%s_LR' % opp[sp]] < 1]

	spots = {}
	for chr, start in zip(d.chr, d.start):
		if chr not in spots:
			spots[chr] = {}
		spots[chr][start] = 1
	
	return spots


def get_ancestral_seq(genomes, chr, spot):
	seq = {}

	for genome in genomes:
		start = spot - 5000
		end = spot + 5000
		seq[genome] = '' 

		seqcall = subprocess.Popen('samtools faidx %s %s:%s-%s' % (genomes[genome], chr, start, end), shell=True, stdout=subprocess.PIPE)
		for l in seqcall.stdout:
			if not re.match('>', l):
				seq[genome] += l.upper().rstrip()

	for genome, sequence in seq.items():
		seq[genome] = list(sequence)

	return seq


def get_lineage_mutations(seq, counts, sp):
	# mutation types
	types = {       'A': {'A': 'AT_AT', 'G': 'AT_GC', 'T': 'AT_AT', 'C': 'AT_GC'},
                        'C': {'A': 'GC_AT', 'G': 'GC_GC', 'T': 'GC_AT', 'C': 'GC_GC'},
                        'T': {'A': 'AT_AT', 'G': 'AT_GC', 'T': 'AT_AT', 'C': 'AT_GC'},
                        'G': {'A': 'GC_AT', 'G': 'GC_GC', 'T': 'GC_AT', 'C': 'GC_GC'},
                }

	for ix, (anc, zf, ltf, dbf) in enumerate(izip(seq['anc'], seq['ZF'], seq['LTF'], seq['DBF'])):
		if sp == 'ZF':
			focal = zf
			others = [ltf, dbf]
		else:
			focal = ltf
			others = [zf, dbf]

		tuple1anc = tuple2anc = tuple1der = tuple2der = 'NN'
                if (ix - 1) >= 0:
                        tuple1anc = seq['anc'][ix - 1] + anc
			tuple1der = seq[sp][ix - 1] + focal
                if (ix + 1) < len(seq['anc']):
                        tuple2anc = anc + seq['anc'][ix + 1]
			tuple2der = focal + seq[sp][ix + 1]

		cpg = False
		for tuple in [tuple1anc, tuple1der, tuple2anc, tuple2der]:
			if tuple in ['CG', 'NG', 'CN']:
				cpg = True
		
		if anc != 'N' and focal in ['A', 'T', 'C', 'G']:
			if not cpg:
				if anc in ['A', 'T']:
					counts[int(ix/100.)]['ancAT'] += 1
				elif anc in ['G', 'C']:
					counts[int(ix/100.)]['ancGC'] += 1
				if focal != anc:
					if focal not in others:
						type = types[anc][focal]
						if type in ['AT_GC', 'GC_AT']:
							counts[int(ix/100.)][type] += 1
	return counts


def main():
	genomes = {	'LTF': '/mnt/gluster/home/sonal.singhal1/reference/LTF_reference.fa',
			'ZF': '/mnt/gluster/home/sonal.singhal1/reference/ZF_reference.fa',
			'DBF': '/mnt/gluster/home/sonal.singhal1/reference/DBF_reference.fa',
			'anc': '/mnt/gluster/home/sonal.singhal1/reference/ancestral_genome.fa'
		}

	sp = 'ZF'
	out = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/bgc_hotspots/unique/unique_%s.hotspots_noCpG.csv' % sp
	o = open(out, 'w')
	o.write('chr,midpoint,species,bin,ancAT,ancGC,AT_GC,GC_AT\n')

	# get unique spots
	spots = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/hotspots.LR_seqldhot.ZF_LTF.csv'
	spots = get_uniq_spots(spots, sp)

	
	for chr in spots:
		for spot in spots[chr]:
			seq = get_ancestral_seq(genomes, chr, spot)
			for sp in ['ZF', 'LTF']:
				counts = {}
				for count in range(0,101):
					counts[count] = {'ancGC': 0, 'ancAT': 0, 'AT_GC': 0, 'GC_AT': 0}

				# identify mutations in focal lineage
				# and get ancestral counts
				counts = get_lineage_mutations(seq, counts, sp)	
				for bin in counts:
					o.write('%s,%s,%s,%s,%s,%s,%s,%s\n' % (chr, spot, sp, bin, counts[bin]['ancAT'], counts[bin]['ancGC'], counts[bin]['AT_GC'], counts[bin]['GC_AT']))
	o.close()

	
if __name__ == "__main__":
    main()
