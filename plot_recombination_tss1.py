import pandas as pd
import re
import numpy as np
from itertools import izip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--chr", help="chromosome for which to run analysis")
args = parser.parse_args()
chr = args.chr

gff = '/mnt/gluster/home/sonal.singhal1/reference/Taeniopygia_guttata.gff'
# long autosomal chromosomes
file = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/without_fam/maps/%s_recombination_bpen100.rm.txt' % chr
out = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/TSS/%s.recombination_tss.csv' % chr
o = open(out, 'w') 
o.write('chr,position,distance_to_TSS,rho\n')

###############

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
	
d = pd.read_csv(file, sep=" ", skiprows=3, header=None, 
	names=['left_snp', 'right_snp', 'meanrho', 'p0.025', 'p0.975'])
rhos = {}
for start, rho in izip(d.left_snp, d.meanrho):
        rhos[start] = rho
max_length = np.max(d.right_snp)
d = ''
starts = sorted(rhos.keys())

def interval_dist(value, start, end):
        middle = int((end - start) / 2.0) + start
        for i in range(start, middle):
                dist = (i - start) * tss[start]
                if (value + 1) < len(starts):
                        if i >= starts[value + 1]:
                                value += 1
                o.write('%s,%s,%s,%s\n' % (chr, i, dist, rhos[starts[value]]))
        for i in range(middle, end):
                dist = (i - end) * tss[end]
                if (value + 1) < len(starts):
                        if i >= starts[value + 1]:
                                value += 1
                o.write('%s,%s,%s,%s\n' % (chr, i, dist, rhos[starts[value]]))
        return value

value = 0
# start of chromosome to the first TSS
for i in range(starts[0], tss_starts[0]):
        dist = (i - tss_starts[0]) * tss[tss_starts[0]]
        if i >= starts[value + 1]:
                value += 1
        o.write('%s,%s,%s,%s\n' % (chr, i, dist, rhos[starts[value]]))
# all the middle intervals
for (start, end) in izip(tss_starts, tss_starts[1:]):
	value = interval_dist(value, start, end)
# from the final TSS to the end of chromosome
for i in range(tss_starts[-1], max_length):
        dist = (i - tss_starts[-1]) * tss[tss_starts[-1]]
        if (value + 1) < len(starts):
                if i >= starts[value + 1]:
                        value += 1
        o.write('%s,%s,%s,%s\n' % (chr, i, dist, rhos[starts[value]]))
o.close()
