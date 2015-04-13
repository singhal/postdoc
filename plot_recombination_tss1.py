import pandas as pd
import re
import numpy as np
from itertools import izip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--chr", help="chromosome for which to run analysis")
parser.add_argument("--sp", help="species for which to run analysis")
args = parser.parse_args()
chr = args.chr
sp = args.sp

gff = '/mnt/gluster/home/sonal.singhal1/reference/Taeniopygia_guttata.gff'
# cpg islands
cpg_file = '/mnt/gluster/home/sonal.singhal1/reference/cpgIslandExt.txt'
# long autosomal chromosomes
file = '/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/maps/%s_recombination_bpen100.txt' % (sp, chr)
out = '/mnt/gluster/home/sonal.singhal1/%s/analysis/TSS/%s.recombination_tss.csv' % (sp, chr)
o = open(out, 'w') 
o.write('chr,position,distance_to_TSS,distance_to_CpG,rho\n')

###############

cpg = []
f = open(cpg_file, 'r')
for l in f:
        d = re.split('\s+', l.rstrip())
        if chr == d[1]:
        	center = int((int(d[3]) + int(d[2])) / 2.0)
        	cpg.append(center)
f.close()

bins = [-10000, -1000, 0, 1000, 10000]

bin_vals = {}
for ix, (i, j) in enumerate(zip(bins, bins[1:])):
	bin_vals[ix+1] = '%s_%s' % (i,j)
bin_vals[0] = '<%s' % bins[0]
bin_vals[len(bins)] = '>%s' % bins[-1]

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

	dists = [x-tss_start for x in cpg]
	abs_dists = [abs(x) for x in dists]
	min_cpg_dist = dists[abs_dists.index(np.min(abs_dists))]

	cpg_bin = np.digitize([min_cpg_dist], bins)[0]
	cpg_bin = bin_vals[cpg_bin]
	print cpg_bin
	tss[ tss_start ]  = {'orientation': orientation, 'cpg_dist' : cpg_bin}
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
                dist = (i - start) * tss[start]['orientation']
                if (value + 1) < len(starts):
                        if i >= starts[value + 1]:
                                value += 1
                o.write('%s,%s,%s,%s,%s\n' % (chr, i, dist, tss[start]['cpg_dist'], rhos[starts[value]]))
        for i in range(middle, end):
                dist = (i - end) * tss[end]['orientation']
                if (value + 1) < len(starts):
                        if i >= starts[value + 1]:
                                value += 1
                o.write('%s,%s,%s,%s,%s\n' % (chr, i, dist, tss[end]['cpg_dist'], rhos[starts[value]]))
        return value

value = 0
# start of chromosome to the first TSS
for i in range(starts[0], tss_starts[0]):
        dist = (i - tss_starts[0]) * tss[tss_starts[0]]['orientation']
        if i >= starts[value + 1]:
                value += 1
        o.write('%s,%s,%s,%s,%s\n' % (chr, i, dist, tss[tss_starts[0]]['cpg_dist'], rhos[starts[value]]))
# all the middle intervals
for (start, end) in izip(tss_starts, tss_starts[1:]):
	value = interval_dist(value, start, end)
# from the final TSS to the end of chromosome
for i in range(tss_starts[-1], max_length):
        dist = (i - tss_starts[-1]) * tss[tss_starts[-1]]['orientation']
        if (value + 1) < len(starts):
                if i >= starts[value + 1]:
                        value += 1
        o.write('%s,%s,%s,%s,%s\n' % (chr, i, dist, tss[tss_starts[-1]]['cpg_dist'], rhos[starts[value]]))
o.close()
