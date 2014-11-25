import pandas as pd
import numpy as np
import sys

ellegren_file = '/mnt/lustre/home/emleffler/finches/zebrafinch/downloads/linkage_map/linkage_map_rec_rate_in_intervals.txt'

longchrs = [ 	'chr1', 'chr1A', 'chr2', 'chr3', 'chr4', 'chr4A', \
				'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', \
				'chr11', 'chr12', 'chr13', 'chr14', 'chr15']

d = pd.read_csv(ellegren_file, sep='\t')
d = d[d.chr.isin(longchrs)]

d['bp_span'] = d.end - d.start
total_bp = d.bp_span.sum()
d = d.sort(['rate'])

half = total_bp / 2.0

sum = 0
for bp, rate in zip(d.bp_span, d.rate):
	sum += bp
	if sum > half:
		print '****'
		print rate
		sys.exit()

