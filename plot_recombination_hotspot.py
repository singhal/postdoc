import re
import pandas as pd
from itertools import izip
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--chr", help="chromosome for which to run analysis")
parser.add_argument("--sp", help="species for which to run analysis")
args = parser.parse_args()
chr = args.chr
sp = args.sp

hot_file = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/spot2kb_flank40kb.seqldhot_validate_hotspots.csv'
rho_file = '/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/maps/%s_recombination_bpen5.txt' % (sp, chr)
# need multiple open files
types = ['zf', 'ltf', 'shared']
out = {}
for type in types:
	out[type] = open('/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/plots/%s_%s.%s_type.csv' % (sp, chr, type), 'w')
	out[type].write('chr,position,dist,rho\n')

# get hotspots in order
d = pd.read_csv(hot_file)
d = d[d.chr == chr]

d['spottype'] = ['NA'] * d.shape[0]
d['shared_start'] = ['NA'] * d.shape[0]
d['shared_length'] = ['NA'] * d.shape[0]
d.ix[ (d.zlk >= 10), 'spottype'] = 'zf'
d.ix[ (d.llk >= 10), 'spottype'] = 'ltf'
d.ix[ (d.zlk >= 10) & (d.llk >= 10), 'spottype'] = 'shared'
if sp == 'ZF':
	d.ix[ (d.spottype == 'shared') | (d.spottype == 'zf'), 'shared_start'] = d.zstart
	d.ix[ (d.spottype == 'shared') | (d.spottype == 'zf'), 'shared_length'] = d.zlength
	d.ix[ (d.spottype == 'ltf'), 'shared_start'] = d.lstart
	d.ix[ (d.spottype == 'ltf'), 'shared_length'] = d.llength
if sp == 'LTF':
	d.ix[ (d.spottype == 'shared') | (d.spottype == 'ltf'), 'shared_start'] = d.lstart
        d.ix[ (d.spottype == 'shared') | (d.spottype == 'ltf'), 'shared_length'] = d.llength
        d.ix[ (d.spottype == 'zf'), 'shared_start'] = d.zstart
        d.ix[ (d.spottype == 'zf'), 'shared_length'] = d.zlength
d = d[d.spottype != 'NA']

d = d.sort(['shared_start'])
# don't want hotspots that are too close to each other
d['dist_diff'] = d.shared_start.diff()
d = d[d.dist_diff > 5000]

rho = pd.read_csv(rho_file, sep=" ", skiprows=3, header=None,
		names=['left_snp', 'right_snp', 'meanrho', 'p0.025', 'p0.975'])

for type, start, length in izip(d.spottype, d.shared_start, d.shared_length):
	midpoint = int(start) + int(int(length) / 2.0)
	block_start = midpoint - 20000
	block_end = midpoint + 20000

	tmp_rho = rho[rho.right_snp >= block_start]
	tmp_rho = tmp_rho[tmp_rho.left_snp <= block_end]

	rhos = {}
	for i, j, meanrho in izip(tmp_rho.left_snp, tmp_rho.right_snp, tmp_rho.meanrho):
		for bp in range(i,j):
			rhos[bp] = meanrho
	
	for bp in range(block_start, block_end):
		if bp in rhos:
			out[type].write('%s,%s,%s,%s\n' % (chr,bp,(bp - midpoint), rhos[bp]))

for type in out:
	out[type].close()
