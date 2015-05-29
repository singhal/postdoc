import pandas as pd
import glob
import re
import numpy as np

file_base = '/mnt/gluster/home/sonal.singhal1/LTF/analysis/LDhelmet/maps/%s.window100000.bpen100.txt'

#longchrs = [    'chr1', 'chr1A', 'chr2', 'chr3', 'chr4', 'chr4A', \
#                                'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', \
#                                'chr11', 'chr12', 'chr13', 'chr14', 'chr15']

chrs = [ 'chr1', 'chr1A', 'chr2', 'chr3',  'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
         'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
         'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28',
         'chrLG2', 'chrLG5', 'chrLGE22', 'chrZ']

rates = []

print 'chr,median_rate,mean_rate'
for chr in chrs:
	if chr not in ['chr16']:
		d = pd.read_csv(file_base % chr)
		print '%s,%s,%s' % (chr, d.rho.median(), d.rho.mean())
		rates += d.rho.tolist()

rates = [x for x in rates if np.isfinite(x)]
print 'all,%s,%s' % (np.median(rates), np.mean(rates))
