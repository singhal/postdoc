import re
import pandas as pd
import os
from itertools import izip
import glob
import numpy as np

files = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/seqldhot_hotspots/*txt')

def parse_file(file):
	file += '.sum'
	zf_lr = None
	if os.path.isfile(file):
		if os.stat(file).st_size > 0:
			d = pd.read_csv(file, sep='\s+', skiprows=1, header=None) 

			d = d[d.X0 > 20000]
			d = d[d.X1 < 30000]
	
			zf_lr = d.X2.max()
	return zf_lr

print 'chr,start,ZF_LR,LTF_LR'
for zf_file in files:
	ltf_file = zf_file.replace('ZF', 'LTF')
	
	zf_lr = parse_file(zf_file)
	ltf_lr = parse_file(ltf_file)

	match = re.search('(chr[\d|Z|A|B]+)_(\d+)', zf_file)
	chr = match.group(1)
	start = match.group(2)

	print '%s,%s,%s,%s' % (chr, start, zf_lr, ltf_lr)
