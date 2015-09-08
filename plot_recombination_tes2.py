import glob
import numpy
import re
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("--sp", help="species for which to run analysis")
args = parser.parse_args()
sp = args.sp

files = glob.glob("/mnt/gluster/home/sonal.singhal1/%s/analysis/TSS/tes/chr*" % sp)
out = '/mnt/gluster/home/sonal.singhal1/%s/analysis/TSS/summary_tes.csv' % (sp)

vals = {}
for file in files:
	f = gzip.open(file,'r')
	header = f.next()
	for l in f:
		d = re.split(',', l.rstrip())

		dist = int(d[3])
		rho = float(d[2])

		if dist not in vals:
			vals[dist] = {'rho': 0, 'num': 0}
		vals[dist]['rho'] += rho
		vals[dist]['num'] += 1

	f.close()

o = open(out, 'w')
o.write('location,average_rho,number_pos\n')
for pos in vals:
	o.write('%s,%s,%s\n' % (pos, vals[pos]['rho'] / float(vals[pos]['num']), vals[pos]['num']))
o.close()
	
