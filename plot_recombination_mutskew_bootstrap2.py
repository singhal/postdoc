import re
import glob
import pandas as pd
import numpy as np
from itertools import izip

dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/mut_skew/'
out = '%ssummary_bootstrap.csv' % (dir)
boot_files = glob.glob('%sboot*csv' % dir)

values = {}
for file in boot_files:
	d = pd.read_csv(file)
	#d = d.rename(columns={'mean': 'mean_rho'})
	for loc, type, mean in izip(d.rho_bin, d.type, d.average_af):
		if loc not in values:
			values[loc] = {}
		if type not in values[loc]:
			values[loc][type] = []
		values[loc][type].append(mean)

o = open(out, 'w')
o.write('rho_bin,type,mean,min,max,p025,p975\n')
for loc in values:
	for type in values[loc]:
		mean = np.mean(values[loc][type])
		min = np.min(values[loc][type])
		max = np.max(values[loc][type])
		p025, p975 = np.percentile(values[loc][type], [2.5, 97.5])

		o.write('%s,%s,%s,%s,%s,%s,%s\n' % (loc, type, mean, min, max, p025, p975))
