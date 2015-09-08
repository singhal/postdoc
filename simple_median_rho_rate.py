import glob
import pandas as pd
import numpy as np
import re

sp = 'ZF'
longchrs = ['chr1', 'chr1A', 'chr2', 'chr3', 'chr4', 'chr4A', 'chr5', 'chr6', \
            'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', \
            'chr15', 'chrZ']
files = glob.glob('/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/maps/*window10000.*' % sp)
longrhos = []
rhos = []

for file in files:
	chr = re.search('(chr[A-Z|0-9]+)', file).group(1)
	if chr in longchrs:
		longrhos += pd.read_csv(file).rho.tolist()
	rhos += pd.read_csv(file).rho.tolist()

rhos = [rho for rho in rhos if np.isfinite(rho)] 
longrhos = [rho for rho in longrhos if np.isfinite(rho)] 

print 'mean (all): %s' % np.mean(rhos)
print 'median (all): %s' % np.median(rhos)
print 'mean (long): %s' % np.mean(longrhos)
print 'median (long): %s' % np.median(longrhos)
