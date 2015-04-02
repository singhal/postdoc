import glob
import numpy
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--sp", help="species for which to run analysis")
args = parser.parse_args()
sp = args.sp

files = glob.glob("/mnt/gluster/home/sonal.singhal1/%s/analysis/TSS/chr*" % sp)
out = '/mnt/gluster/home/sonal.singhal1/%s/analysis/TSS/cpg/summary_cpg.csv' % sp

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
				cpg[dist][cpg_bin] = {'rho': 0, 'num': 0}
			cpg[dist][cpg_bin]['rho'] += rho
			cpg[dist][cpg_bin]['num'] += 1
	f.close()

o = open(out, 'w')
o.write('location,cpg_bin,average_rho,number_pos\n')
for pos in cpg:
	for cpg_bin in cpg[pos]:
		o.write('%s,%s,%s,%s\n' % (pos, cpg_bin, cpg[pos][cpg_bin]['rho'] / float(cpg[pos][cpg_bin]['num']), cpg[pos][cpg_bin]['num']))
o.close()
	
