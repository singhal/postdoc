import pandas as pd
import numpy as np
import random
import scipy as sci
import re
import glob
import gzip

files = glob.glob("/mnt/gluster/home/sonal.singhal1/ZF/analysis/TSS/chr*")
dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/TSS/bootstrap/'

values = {}
for file in files:
        f = gzip.open(file,'r')
        header = f.next()
        for l in f:
                d = re.split(',', l.rstrip())

                tss_bin = ((int(d[3]) / 500) + 1) * 500
		cpg_bin = ((int(d[4]) / 500) + 1) * 500
		rho = float(d[2])

                if abs(tss_bin) <= 1e5 and abs(cpg_bin) <= 1e5:
			if tss_bin not in values:
				values[tss_bin] = {}
			if cpg_bin not in values[tss_bin]:
				values[tss_bin][cpg_bin] = []
			values[tss_bin][cpg_bin].append(rho)
	f.close()

for bootstrap in range(0,100):
	out = '%sbootstrap%s.csv' % (dir, bootstrap)
	o = open(out, 'w')
	o.write('tss_dist,cpg_dist,mean,num\n')
	for tss_dist in values:
		for cpg_dist in values[tss_dist]:
			original = np.array(values[tss_dist][cpg_dist])
			resample = np.floor(np.random.rand(len(original))*len(original)).astype(int)
			resampled_mean = np.mean(original[resample])
			o.write('%s,%s,%s,%s\n' % (tss_dist, cpg_dist, resampled_mean, len(original)))
	o.close()
