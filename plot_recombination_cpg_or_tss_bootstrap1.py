import glob
import numpy as np
import re
import argparse
import gzip
import random
import scipy as sci

parser = argparse.ArgumentParser()
parser.add_argument("--sp", help="species for which to run analysis")
parser.add_argument("--type", help="tss or cpg?")
args = parser.parse_args()
tss_cpg = args.type
sp = args.sp

files = glob.glob("/mnt/gluster/home/sonal.singhal1/%s/analysis/TSS/chr*" % sp)
outdir = '/mnt/gluster/home/sonal.singhal1/%s/analysis/TSS/bootstrap_%s/' % (sp, tss_cpg)

if tss_cpg == 'cpg':
	pos = 4
	binnum = 6
	binname = 'tss_bin'
elif tss_cpg == 'tss':
	pos = 3
	binnum = 5
	binname = 'cpg_bin'

vals = {}
for file in files:
	f = gzip.open(file,'r')
	header = f.next()
	for l in f:
		d = re.split(',', l.rstrip())

		dist = int(d[pos])
		bin = d[binnum]
		rho = float(d[2])

		if abs(dist) <= 1e5:
			if dist not in vals:
				vals[dist] = {}
			if bin not in vals[dist]:
				vals[dist][bin] = []
			vals[dist][bin].append(rho)
	f.close()


for bootstrap in range(0,100):
        out = '%sbootstrap%s.csv' % (outdir, bootstrap)
        o = open(out, 'w')
        o.write('location,%s,mean,num\n' % binname)
        for dist in vals:
                for bin in vals[dist]:
                        original = np.array(vals[dist][bin])
                        resample = np.floor(np.random.rand(len(original))*len(original)).astype(int)
                        resampled_mean = np.mean(original[resample])
                        o.write('%s,%s,%s,%s\n' % (dist, bin, resampled_mean, len(original)))
        o.close()
