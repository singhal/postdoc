import glob
import numpy
import re
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("--sp", help="species for which to run analysis")
args = parser.parse_args()
sp = args.sp

files = glob.glob("/mnt/gluster/home/sonal.singhal1/%s/analysis/TSS/chr*" % sp)
out = '/mnt/gluster/home/sonal.singhal1/%s/analysis/TSS/summary_cpg_dist_to_TSS.csv' % sp

cpg = {}
for file in files:
	f = gzip.open(file,'r')
	header = f.next()
	for l in f:
		d = re.split(',', l.rstrip())

		dist = int(d[4])
		bin = d[6]
		rho = float(d[2])

		if abs(dist) <= 1e5:
			if dist not in cpg:
				cpg[dist] = {}
			if bin not in cpg[dist]:
				cpg[dist][bin] = {'rho': 0, 'num': 0}
			cpg[dist][bin]['rho'] += rho
			cpg[dist][bin]['num'] += 1
	f.close()

o = open(out, 'w')
o.write('location,tss_bin,average_rho,number_pos\n')
for pos in cpg:
	for bin in cpg[pos]:
		o.write('%s,%s,%s,%s\n' % (pos, bin, cpg[pos][bin]['rho'] / float(cpg[pos][bin]['num']), cpg[pos][bin]['num']))
o.close()
	
