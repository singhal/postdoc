import glob
import numpy
import re
import argparse
import random
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--sp", help="species for which to run analysis")
args = parser.parse_args()
sp = args.sp

files = glob.glob("/mnt/gluster/home/sonal.singhal1/%s/analysis/TSS/chr*" % sp)

cpg = {}
for file in files:
	f = open(file,'r')
	header = f.next()
	for l in f:
		d = re.split(',', l.rstrip())

		dist = int(d[2])
		cpg_bin = d[3]
		rho = float(d[4])

		if abs(dist) <= 1e5:
			if dist not in cpg:
				cpg[dist] = {}
			if cpg_bin not in cpg[dist]:
				cpg[dist][cpg_bin] = []
			cpg[dist][cpg_bin].append(rho)
	f.close()

for bootstrap in range(0,100):
	out = '/mnt/gluster/home/sonal.singhal1/%s/analysis/TSS/cpg/bootstrap%s.csv' % (sp, bootstrap)
	o = open(out, 'w')
	o.write('location,cpg_bin,mean\n')
	for pos in cpg:
		for cpg_bin in cpg[pos]:
			original = np.array(cpg[pos][cpg_bin])
                	resample = np.floor(np.random.rand(len(original))*len(original)).astype(int)
                	resampled_mean = np.mean(original[resample])
			o.write('%s,%s,%s\n' % (pos, cpg_bin, resampled_mean))
	o.close()
	
