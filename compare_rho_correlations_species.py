import pandas as pd
from scipy.stats.stats import pearsonr
import numpy as np

longchrs = [	'chr1', 'chr1A', 'chr2', 'chr3', 'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', \
                'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15']
# no chr16 because failed for LTF
chrs = [ 'chr1', 'chr1A', 'chr1B', 'chr2', 'chr3',  'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
         'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr17', 'chr18',
         'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28',
         'chrLG2', 'chrLG5', 'chrLGE22', 'chrZ' ]

out = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/compare_rho_two_species.csv'
o = open(out, 'w')
o.write('chr,window,corr,pval,deviation\n')
for chr in longchrs:
	for window in [10000, 100000, 1000000, 5000000]:
		zf = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/maps/%s.window%s.bpen100.txt' % (chr, window)
		ltf = zf.replace('ZF', 'LTF')
	
		zfd = pd.read_csv(zf)
		ltfd = pd.read_csv(ltf)

		rho1 = []
		rho2 = []
		dev = []
		for a, b in zip(zfd.rate, ltfd.rate):
			if np.isfinite(a) and np.isfinite(b):
				rho1.append(a)
				rho2.append(b)
				deviation = abs(1 - (a/b))
				if np.isfinite(deviation):
					dev.append( deviation )			

		if len(rho1) > 2:
			r, pval = pearsonr(rho1, rho2)
			o.write('%s,%s,%.3f,%.3f,%.3f\n' % (chr, window, r, pval, np.mean(dev)))
o.close()
