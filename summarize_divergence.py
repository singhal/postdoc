import re
import glob
import pandas as pd
import numpy as np

types = ['ZF_LTF', 'ZF_LTFh', 'ZF_LTFa', 'LTFa_LTFh', 'ZF_DBF', 'LTF_DBF']

for type in types:
	files = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/analysis/pop_gen/*%s.divergence*' % type)
	files = [x for x in files if not re.search('chrZ', x)]

	num_sites = 0
	dxy_sum = 0

	for file in files:
		d = pd.read_csv(file)

		num_sites += d.num_sites.sum()
		dxy_sum += (d.num_sites * d.dxy).sum()

	dxy = dxy_sum / float(num_sites)

	print '%s %s %.4f' % (type, num_sites, dxy)
