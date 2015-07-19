import glob
import pandas as pd
import numpy as np

sp = 'LTF'
files = glob.glob('/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/maps/*window10000.*' % sp)
rhos = []

for file in files:
	rhos += pd.read_csv(file).rho.tolist()

rhos = [rho for rho in rhos if np.isfinite(rho)] 
#print rhos
print np.median(rhos)
