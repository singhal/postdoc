import glob
import numpy as np
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--sp", help="species for which to run analysis")
parser.add_argument("--type", help="tss or cpg?")
args = parser.parse_args()
tss_cpg = args.type
sp = args.sp

dir = '/mnt/gluster/home/sonal.singhal1/%s/analysis/TSS/bootstrap_%s/' % (sp, tss_cpg)
out = '/mnt/gluster/home/sonal.singhal1/%s/analysis/TSS/summary_%s_bootstrap.csv' % (sp, tss_cpg)

files = glob.glob('%sboot*csv' % dir)

vals = {}
for file in files:
	f = open(file,'r')
	header = f.next()
	for l in f:
		d = re.split(',', l.rstrip())
		if d[0] not in vals:
			vals[d[0]] = {}
		if d[1] not in vals[d[0]]:
			vals[d[0]][d[1]] = {'num': d[3], 'vals': []}
		vals[d[0]][d[1]]['vals'].append(float(d[2]))
	f.close()


o = open(out, 'w')
o.write('location,bin,num,mean,min,max,p025,p975\n')
for loc in vals:
        for bin in vals[loc]:
                mean = np.mean(vals[loc][bin]['vals'])
                min = np.min(vals[loc][bin]['vals'])
                max = np.max(vals[loc][bin]['vals'])
                p025, p975 = np.percentile(vals[loc][bin]['vals'], [2.5, 97.5])

                o.write('%s,%s,%s,%s,%s,%s,%s,%s\n' % (loc, bin, vals[loc][bin]['num'], mean, min, max, p025, p975))
o.close()
