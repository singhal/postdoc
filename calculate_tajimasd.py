import re
import pandas as pd
from itertools import izip
import math

sp = 'ZF'

chr_lengths = { 'chr10': 20806668, 'chr11': 21403021, 'chr12': 21576510, 'chr13': 16962381,
                'chr14': 16419078, 'chr15': 14428146, 'chr16': 9909, 'chr17': 11648728,
                'chr18': 11201131, 'chr19': 11587733, 'chr1A': 73657157, 'chr1B': 1083483,
                'chr1': 118548696, 'chr20': 15652063, 'chr21': 5979137, 'chr22': 3370227,
                'chr23': 6196912, 'chr24': 8021379, 'chr25': 1275379, 'chr26': 4907541,
                'chr27': 4618897, 'chr28': 4963201, 'chr2': 156412533, 'chr3': 112617285,
                'chr4A': 20704505, 'chr4': 69780378, 'chr5': 62374962, 'chr6': 36305782,
                'chr7': 39844632, 'chr8': 27993427, 'chr9': 27241186, 'chrLG2': 109741,
                'chrLG5': 16416, 'chrLGE22': 883365, 'chrZ': 72861351}

if sp == 'ZF':
        ninds = {}
        for chr in chr_lengths:
                ninds[chr] = 38
        ninds['chrZ'] = 28
if sp == 'LTF':
        ninds = {}
        for chr in chr_lengths:
                ninds[chr] = 40
        ninds['chrZ'] = 32

pifile = '/mnt/gluster/home/sonal.singhal1/%s/analysis/pop_gen/pi.csv' % sp
thetafile = '/mnt/gluster/home/sonal.singhal1/%s/analysis/pop_gen/wattersons_theta.csv' % sp
outfile = '/mnt/gluster/home/sonal.singhal1/%s/analysis/pop_gen/tajimasd.csv' % sp

pi = pd.read_csv(pifile)
theta = pd.read_csv(thetafile)

d = pd.merge(pi, theta, on=['chr', 'index', 'start'])

tajd = []
for chr, pi, theta, seqlength in izip(d.chr, d.pi, d.watterson_theta, d.seq_length_x):
	if not pd.isnull(pi) and not pd.isnull(theta):
		a1 = 0
		a2 = 0

		nind = ninds[chr]	
	
		for i in range(1, nind):
			a1 += 1 / float(i)
			a2 += 1 / float(i ** 2)
		S = theta * a1 * seqlength
		b1 = (nind + 1) / (3.0 * (nind - 1))
		b2 = 2 * (nind ** 2 + nind + 3) / (9.0 * nind * (nind - 1))
		c1 = b1 - 1 / float(a1)
		c2 = b2 - (nind + 2) / (float(a1) * nind) + a2 / float(a1 ** 2)
		e1 = c1 / float(a1)
		e2 = c2 / float(a1 ** 2 + a2)
		variance = math.sqrt(e1 * S + e2 * S * (S - 1))
		pilength = pi * seqlength
		tajd_tmp = (pilength - S / float(a1) ) / float(variance)
		tajd.append(tajd_tmp) 
	else:
		tajd.append('NA')
d['tajimasd'] = tajd
d.to_csv(outfile, index=False)
