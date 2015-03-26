import glob
import numpy
import re

files = glob.glob("/mnt/gluster/home/sonal.singhal1/LTF/analysis/TSS/chr*")
out = '/mnt/gluster/home/sonal.singhal1/LTF/analysis/TSS/summary_tss.csv'

tss = {}
for file in files:
	f = open(file,'r')
	header = f.next()
	for l in f:
		d = re.split(',', l.rstrip())

		d[2] = int(d[2])
		d[3] = float(d[3])

		if abs(d[2]) <= 1e6:
			if d[2] not in tss:
				tss[d[2]] = {'rho': 0, 'num': 0}
			tss[d[2]]['rho'] += d[3]
			tss[d[2]]['num'] += 1

o = open(out, 'w')
o.write('location,average_rho,number_pos\n')
for pos in tss:
	o.write('%s,%s,%s\n' % (pos, tss[pos]['rho'] / float(tss[pos]['num']), tss[pos]['num']))
o.close()
	
