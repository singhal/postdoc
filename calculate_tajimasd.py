import re
import pandas as pd
from itertools import izip
import math

species = 'ZF'
nind = 38

pifile = '/mnt/gluster/home/sonal.singhal1/%s/analysis/pop_gen/pi.csv' % species
thetafile = '/mnt/gluster/home/sonal.singhal1/%s/analysis/pop_gen/wattersons_theta.csv' % species
outfile = '/mnt/gluster/home/sonal.singhal1/%s/analysis/pop_gen/tajimasd.csv' % species

pi = pd.read_csv(pifile)
theta = pd.read_csv(thetafile)

d = pd.merge(pi, theta, on=['chr', 'index', 'start', 'end'])

tajd = []
for pi, theta, seqlength in izip(d.pi, d.watterson_theta, d.seq_length):
	if not pd.isnull(pi) and not pd.isnull(theta):
		a1 = 0
		a2 = 0
	
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
