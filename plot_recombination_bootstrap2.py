import re
import glob
import pandas as pd
import numpy as np

dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/TSS/bootstrap/'
out = '%ssummary_bootstrap.csv' % (dir)
boot_files = glob.glob('%sboot*csv' % dir)

values = {}
for file in boot_files:
	print file
	d = pd.read_csv(file)
	d = d.rename(columns={'mean': 'mean_rho'})
	for tss, cpg, mean, num in zip(d.tss_dist, d.cpg_dist, d.mean_rho, d.num):
		if tss not in values:
			values[tss] = {}
		if cpg not in values[tss]:
			values[tss][cpg] = {'vals': [], 'num': num}
		values[tss][cpg]['vals'].append(mean)

o = open(out, 'w')
o.write('tss_dist,cpg_dist,num,mean,min,max,p025,p975\n')
for tss in values:
	for cpg in values[tss]:
		mean = np.mean(values[tss][cpg]['vals'])
		min = np.min(values[tss][cpg]['vals'])
		max = np.max(values[tss][cpg]['vals'])
		p025, p975 = np.percentile(values[tss][cpg]['vals'], [2.5, 97.5])

		o.write('%s,%s,%s,%s,%s,%s,%s,%s\n' % (tss, cpg, values[tss][cpg]['num'], mean, min, max, p025, p975))
