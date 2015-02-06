import re
import pandas as pd
import numpy as np

chrs = [ 'chr1', 'chr1A', 'chr1B', 'chr2', 'chr3',  'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
         'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
         'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28',
         'chrLG2', 'chrLG5', 'chrLGE22']

def get_average(file):
	d = pd.read_csv(file)
	avg = (d.breaks - d.big_breaks ) / [float(x) for x in d.het_sites]
	return np.median( avg )


avg_change = []
for chr in chrs:
	pre = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/breaks/recombinationbreaks.%s.hapi.csv' % chr
	post = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/breaks/recombinationbreaks.%s.hapi.noswitch.csv' % chr

	pre = get_average(pre)
	post = get_average(post)
	change = (post - pre) / pre 
	print '%s %.4f %.4f %.1f' % (chr, pre, post, change)
	avg_change.append(change)

print 'Average change is : %.3f' % (np.mean(filter(lambda x: np.isfinite(x), avg_change)))
