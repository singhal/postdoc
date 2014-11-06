import re
from scipy import stats
import numpy as np
import pandas as pd

d_file = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/simulations/compare_sim_to_real1.csv'
d = pd.read_csv(d_file)

print 'bpen,rho,r,median_diff,residuals'
groups = d.groupby(d.bpen)
for bpen, group in groups:
	grouped = group.groupby(d.rho_strength)
	for rho, group2 in grouped:
		r, pval = stats.pearsonr(group2.actual_rec, group2.estimated_rec)
		
		residuals = sum(group2.actual_rec - group2.estimated_rec)**2
		print '%s,%s,%.2f,%.2f,%.2f' % (bpen, rho, r, np.median(group2.diffratio), residuals)
