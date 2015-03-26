import glob
import numpy
import re

files = glob.glob("/mnt/gluster/home/sonal.singhal1/ZF/analysis/cpg/chr*")
out = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/cpg/summary_cpg.csv'

cpg = {}
for file in files:
	f = open(file,'r')
	header = f.next()
	for l in f:
		d = re.split(',', l.rstrip())

		d[2] = int(d[2])
		d[3] = float(d[3])

		if abs(d[2]) <= 1e6:
			if d[2] not in cpg:
				cpg[d[2]] = {'rho': 0, 'num': 0}
			cpg[d[2]]['rho'] += d[3]
			cpg[d[2]]['num'] += 1
	f.close()

o = open(out, 'w')
o.write('location,average_rho,number_pos\n')
for pos in cpg:
	o.write('%s,%s,%s\n' % (pos, cpg[pos]['rho'] / float(cpg[pos]['num']), cpg[pos]['num']))
o.close()
	
