import pandas as pd
import glob
import re
import numpy as np

file_base = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/maps/%s.window100000.bpen100.txt'

longchrs = [    'chr1', 'chr1A', 'chr2', 'chr3', 'chr4', 'chr4A', \
                                'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', \
                                'chr11', 'chr12', 'chr13', 'chr14', 'chr15']

rates = []

print 'chr,median_rate,mean_rate'
for chr in longchrs:
	d = pd.read_csv(file_base % chr)
	print '%s,%s,%s' % (chr, d.rho.median(), d.rho.mean())
	rates += d.rho.tolist()

rates = [x for x in rates if np.isfinite(x)]
print 'all,%s,%s' % (np.median(rates), np.mean(rates))
