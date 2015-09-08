import glob
import numpy
import re
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("--sp", help="species for which to run analysis")
parser.add_argument("--type", help="tss or cpg?")
args = parser.parse_args()
tss_cpg = args.type
sp = args.sp

files = glob.glob("/mnt/gluster/home/sonal.singhal1/%s/analysis/TSS/chr*" % sp)
out = '/mnt/gluster/home/sonal.singhal1/%s/analysis/TSS/summary_%s.csv' % (sp, tss_cpg)

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
				vals[dist][bin] = {'rho': 0, 'num': 0}
			vals[dist][bin]['rho'] += rho
			vals[dist][bin]['num'] += 1
	f.close()

o = open(out, 'w')
o.write('location,%s,average_rho,number_pos\n' % binname)
for pos in vals:
	for bin in vals[pos]:
		o.write('%s,%s,%s,%s\n' % (pos, bin, vals[pos][bin]['rho'] / float(vals[pos][bin]['num']), vals[pos][bin]['num']))
o.close()
	
