import re
import pandas as pd
from itertools import izip
import numpy as np
import argparse

# add tail end

parser = argparse.ArgumentParser()
parser.add_argument("--chr", help="chromosome for which to run analysis")
args = parser.parse_args()
chr = args.chr

cpg_file = '/mnt/gluster/home/sonal.singhal1/reference/cpgIslandExt.txt'
rho_file = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/maps/%s_recombination_bpen100.txt' % chr
out = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/cpg/%s.cpg_recombination.csv' % chr
o = open(out, 'w')
o.write('chr,position,distance,rho\n')

cpg = []
f = open(cpg_file, 'r')
for l in f:
	d = re.split('\s+', l.rstrip())
	if d[1] == chr:
		cpg_length = int(d[3]) - int(d[2])
		# getting rid of very short cpgs
		if cpg_length >= 300:
			center = int((int(d[3]) + int(d[2])) / 2.0)
			cpg.append(center)
f.close()

rhos = {}
d = pd.read_csv(rho_file, sep=" ", skiprows=3, header=None,
		names=['left_snp', 'right_snp', 'meanrho', 'p0.025', 'p0.975'])
for start, rho in izip(d.left_snp, d.meanrho):
	rhos[start] = rho
max_length = np.max(d.right_snp)
d = ''
starts = sorted(rhos.keys())

def interval_dist(value, start, end):
	middle = int((end - start) / 2.0) + start
	for i in range(start, middle):
		dist = i - start
		if (value + 1) < len(starts):
			if i >= starts[value + 1]:
                		value += 1
        	o.write('%s,%s,%s,%s\n' % (chr, i, dist, rhos[starts[value]]))
	for i in range(middle, end):
		dist = i - end
		if (value + 1) < len(starts):
			if i >= starts[value + 1]:
                		value += 1
        	o.write('%s,%s,%s,%s\n' % (chr, i, dist, rhos[starts[value]]))
	return value

value = 0
for i in range(starts[0], cpg[0]):
	dist = i - cpg[0]
	if i >= starts[value + 1]:
		value += 1
	o.write('%s,%s,%s,%s\n' % (chr, i, dist, rhos[starts[value]]))
for (start, end) in zip(cpg, cpg[1:]):
	value = interval_dist(value, start, end)
for i in range(cpg[-1], max_length):
	dist = i - cpg[-1]
	if (value + 1) < len(starts):
	        if i >= starts[value + 1]:
        	        value += 1
        o.write('%s,%s,%s,%s\n' % (chr, i, dist, rhos[starts[value]]))
o.close()
