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
out = '/mnt/gluster/home/sonal.singhal1/%s/analysis/TSS/%s.recombination_tss_cpg.csv' % (sp, chr)
o = open(out, 'w') 

###############

cpg_starts = []
f = open(cpg_file, 'r')
for l in f:
        d = re.split('\s+', l.rstrip())
        if chr == d[1]:
        	center = int((int(d[3]) + int(d[2])) / 2.0)
        	cpg_starts.append(center)
f.close()
cpg_starts = sorted(cpg_starts)


# CpG bins
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
	
	dists = [x-tss_start for x in cpg_starts]
        abs_dists = [abs(x) for x in dists]
        min_cpg_dist = dists[abs_dists.index(np.min(abs_dists))]

        cpg_bin = np.digitize([min_cpg_dist], bins)[0]
        cpg_bin = bin_vals[cpg_bin]

	tss[ tss_start ]  = {'orientation': orientation, 'cpg_dist' : cpg_bin}
tss_starts = sorted(tss.keys())

d = pd.read_csv(file, sep=" ", skiprows=3, header=None, 
	names=['left_snp', 'right_snp', 'meanrho', 'p0.025', 'p0.975'])
rhos = {}
for start, rho in izip(d.left_snp, d.meanrho):
        rhos[start] = rho
max_length = np.max(d.right_snp)
d = ''
rho_starts = sorted(rhos.keys())

tss_value = 0
cpg_value = 0
rho_value = 0

o.write('chr,position,rho,TSS_dist,CpG_dist,CpG_dist_TSS\n')
for bp in range(rho_starts[0], max_length):
	if (rho_value + 1) < len(rho_starts):
		if bp > rho_starts[rho_value + 1]:
			rho_value += 1
	if (cpg_value + 1) < len(cpg_starts):
		# switch to next island because the position is greater than the midpoint between
		#	two adjacent islands
		if bp > int(0.5 * (cpg_starts[cpg_value] + cpg_starts[cpg_value + 1])):
			cpg_value += 1
	if (tss_value + 1) < len(tss_starts):
		# switch to the next TSS because the pos is greater than the midpoint between
		#	two adjacent tss
		if bp > int(0.5 * (tss_starts[tss_value] + tss_starts[tss_value + 1])):
        	        tss_value += 1

	rho_val = rhos[rho_starts[rho_value]]
	tss_dist = (bp - tss_starts[tss_value]) * tss[tss_starts[tss_value]]['orientation']
	cpg_dist = abs(bp - cpg_starts[cpg_value]) 
	tss_cpg = tss[tss_starts[tss_value]]['cpg_dist']

	o.write('%s,%s,%s,%s,%s,%s\n' % (chr, bp, rho_val, tss_dist, cpg_dist, tss_cpg))
o.close()
