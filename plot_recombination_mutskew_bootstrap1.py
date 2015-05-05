import pandas as pd
import glob
import numpy as np
import re
import argparse
import random

parser = argparse.ArgumentParser()
parser.add_argument("--sp", help="species for which to run analysis")
args = parser.parse_args()
sp = args.sp

rhos = []
chrs = ['chr1', 'chr1A', 'chr2', 'chr3', 'chr4', 'chr4A', 'chr5', 'chr6', 'chr8', \
	'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15']
for chr in chrs:
	d = pd.read_csv('/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/maps/%s.window10000.bpen100.txt' % (sp, chr))
	rhos += d.rho.tolist()
rhos = [x for x in rhos if np.isfinite(x)]

bins = np.percentile(rhos, range(0,101))

files = glob.glob("/mnt/gluster/home/sonal.singhal1/%s/analysis/mut_skew/chr*" % sp)

types = ['AT_GC', 'GC_AT', 'GC_GC', 'AT_AT']
skew = {}
for file in files:
	f = open(file,'r')
	header = f.next()
	for l in f:
		d = re.split(',', l.rstrip())

		rho = float(d[5])
		af = float(d[3])
		type = d[2]
		rho_bin = np.digitize([rho], bins)[0]

		if af < 1 and af > 0:
			if rho_bin not in skew:
				skew[rho_bin] = {}
				for t in types:
					skew[rho_bin][t] = {'rho': 0, 'num': 0, 'af': []}
			skew[rho_bin][type]['rho'] += rho
			skew[rho_bin][type]['num'] += 1
			skew[rho_bin][type]['af'].append(af)
	f.close()

rho_vals = {}
for ix, (lo, hi) in enumerate(zip(bins, bins[1:])):
	bin_ix = ix + 1
	rho_vals[bin_ix] = '%.4f_%.4f' % (lo, hi)
rho_vals[0] = '<%.4f' % (bins[0])
rho_vals[len(bins)] = '>%.4f' % (bins[-1])

for bootstrap in range(0,100):
	out = '/mnt/gluster/home/sonal.singhal1/%s/analysis/mut_skew/bootstrap%s.csv' % (sp, bootstrap)
	o = open(out, 'w')
	o.write('rho_bin,rho_val,type,average_af,average_rho,number_pos\n')
	for rho_bin in skew:
		for type in skew[rho_bin]:
			if skew[rho_bin][type]['num'] > 0:
				original = np.array(skew[rho_bin][type]['af'])
                		resample = np.floor(np.random.rand(len(original))*len(original)).astype(int)
                		resampled_mean = np.mean(original[resample])

				avg_rho = skew[rho_bin][type]['rho'] / float( skew[rho_bin][type]['num'] )
				o.write('%s,%s,%s,%s,%s,%s\n' % (rho_bin, rho_vals[rho_bin], type, resampled_mean, avg_rho, skew[rho_bin][type]['num'] ))
o.close()
