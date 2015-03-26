import re
import glob
import pandas as pd
import numpy as np

dir = '/mnt/gluster/home/sonal.singhal1/LTF/analysis/TSS/'
out = '%ssummary_bootstrap.csv' % (dir)
boot_files = glob.glob('%sboot*csv' % dir)

values = {}
for file in boot_files:
	d = pd.read_csv(file)
	d = d.rename(columns={'mean': 'mean_rho'})
	for loc, mean in zip(d.location, d.mean_rho):
		if loc not in values:
			values[loc] = []
		values[loc].append(mean)

o = open(out, 'w')
o.write('loc,mean,min,max,p025,p975\n')
for loc in values:
	mean = np.mean(values[loc])
	min = np.min(values[loc])
	max = np.max(values[loc])
	p025, p975 = np.percentile(values[loc], [2.5, 97.5])

	o.write('%s,%s,%s,%s,%s,%s\n' % (loc, mean, min, max, p025, p975))
