import pandas as pd
from itertools import izip

rho_min = 0.001
rho_max = 0.6
        
chrs = ['chr1', 'chr1A', 'chr2', 'chr3', 'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
        'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chrZ']

match = {}
for sp in ['ZF', 'LTF']:
	if sp == 'ZF':
		rho_dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/without_fam/maps/'
	if sp == 'LTF':
		rho_dir = '/mnt/gluster/home/sonal.singhal1/LTF/analysis/LDhelmet/old/maps/'
	
	seq = 0	

	for chr in chrs:
		file = '%s/%s.window100000.bpen100.rm.txt' % (rho_dir, chr)
		d = pd.read_csv(file)

		d = d.rename(columns={'rate': 'rho_rate'})
		d = d[d.rho_rate >= rho_min]
		d = d[d.rho_rate <= rho_max]

		for chr, start, end in izip(d.chr, d.window_start, d.window_end):
			if chr not in match:
				match[chr] = {}
			if start not in match[chr]:
				match[chr][start] = {end: 0}
			match[chr][start][end] += 1
			seq += (end - start)
	print '%s %s' % (sp, seq)

shared = 0
for chr in match:
	for start in match[chr]:
		for end in match[chr][start]:
			if match[chr][start][end] == 2:
				shared += (end - start)
print '%s %s' % ('shared', shared)	
	

