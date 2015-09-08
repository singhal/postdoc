import pandas as pd
import re
import numpy as np
from itertools import izip
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("--chr", help="chromosome for which to run analysis")
parser.add_argument("--sp", help="species for which to run analysis")
args = parser.parse_args()
chr = args.chr
sp = args.sp

gff = '/mnt/gluster/home/sonal.singhal1/reference/Taeniopygia_guttata.gff'
# long autosomal chromosomes
file = '/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/maps/%s_recombination_bpen100.txt' % (sp, chr)
out = '/mnt/gluster/home/sonal.singhal1/%s/analysis/TSS/tes/%s.recombination_tes.csv.gz' % (sp, chr)
o = gzip.open(out, 'w') 

###############

d = pd.read_csv(gff, sep='\t', header=0, names=['chr', 'type', 'cds_mrna', 
						'start', 'stop', 'score', 
						'orientation', 'codon_pos', 'id'])
chr_short = chr.replace('chr','')
d = d[d.chr == chr_short]
gdata = d[ d.cds_mrna == 'mRNA' ].set_index('id').to_dict()
genes = d[ d.cds_mrna == 'mRNA' ].id.unique()

tes = {}
for gene in genes:
	orientation = 1
	if gdata['orientation'][gene] == '+':
		tes_end = gdata['stop'][gene]
	elif gdata['orientation'][gene]  == '-':
		tes_end = gdata['start'][gene]
		orientation = -1

	tes[ tes_end ]  = orientation
tes_ends = sorted(tes.keys())

d = pd.read_csv(file, sep=" ", skiprows=3, header=None, 
	names=['left_snp', 'right_snp', 'meanrho', 'p0.025', 'p0.975'])
rhos = {}
for start, rho in izip(d.left_snp, d.meanrho):
        rhos[start] = rho
max_length = np.max(d.right_snp)
d = ''
rho_starts = sorted(rhos.keys())

tes_value = 0
rho_value = 0

o.write('chr,position,rho,TES_dist\n')
for bp in range(rho_starts[0], max_length):
	if (rho_value + 1) < len(rho_starts):
		if bp > rho_starts[rho_value + 1]:
			rho_value += 1
	if (tes_value + 1) < len(tes_ends):
		# switch to the next TSS because the pos is greater than the midpoint between
		#	two adjacent tss
		if bp > int(0.5 * (tes_ends[tes_value] + tes_ends[tes_value + 1])):
        	        tes_value += 1

	rho_val = rhos[rho_starts[rho_value]]
	tes_dist = (bp - tes_ends[tes_value]) * tes[tes_ends[tes_value]]

	o.write('%s,%s,%s,%s\n' % (chr, bp, rho_val, tes_dist))
o.close()
print 'done!'
