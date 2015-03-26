import pandas as pd
import numpy as np
import random
import scipy as sci
import re
import glob

files = glob.glob("/mnt/gluster/home/sonal.singhal1/LTF/analysis/TSS/chr*")
dir = '/mnt/gluster/home/sonal.singhal1/LTF/analysis/TSS/'

values = {}
for file in files:
        f = open(file,'r')
        header = f.next()
        for l in f:
                d = re.split(',', l.rstrip())

                d[2] = int(d[2])
                d[3] = float(d[3])

                if abs(d[2]) <= 1e5:
                        if d[2] not in values:
                                values[d[2]] = []
                        values[d[2]].append(d[3])

for bootstrap in range(0,100):
	out = '%sbootstrap%s.csv' % (dir, bootstrap)
	o = open(out, 'w')
	o.write('location,mean\n')
	for loc in values:
		original = np.array(values[loc])
		resample = np.floor(np.random.rand(len(original))*len(original)).astype(int)
		resampled_mean = np.mean(original[resample])
		o.write('%s,%s\n' % (loc, resampled_mean))
	o.close()
