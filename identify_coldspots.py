import re
import pandas as pd
import glob
import numpy as np

files = glob.glob('/mnt/gluster/home/sonal.singhal1/LTF/analysis/LDhelmet/putative_hotspots/*block2000_flank40000*')
out = '/mnt/gluster/home/sonal.singhal1/for_ellen/ldhelmet_unvalidated_coldspots.LTF.csv'
o = open(out, 'w')
o.write('chr,spot_start,spot_end,flank_rate,spot_heat\n')

# borrowed from 
# http://nbviewer.ipython.org/github/goldenhelix/cftr-variant-classification-analysis/blob/master/CFTR%20Data%20Munging.ipynb
def find_intervals(intervals):
    cur_start, cur_stop = intervals[0]
    for next_start, next_stop in intervals[1:]:
        if cur_stop <= next_start:
            yield cur_start, cur_stop
            cur_start, cur_stop = next_start, next_stop
        else:
            cur_stop = next_stop
    yield cur_start, cur_stop

for file in files:
	d = pd.read_csv(file)
	d = d[d.rate_ratio < 0.1]
	d = d[d.flank_rate > 0.001]
	d = d[d.flank_rate < 0.6]
	d = d[np.isfinite(d.block_rate)]
	
	# for chr, block_start, flank_rate, spot_heat in zip(d.chr, d.block_start, d.flank_rate, d.rate_ratio):
	#	o.write('%s,%s,%s,%s\n' % (chr, block_start, flank_rate, spot_heat))
	
	intervals = []
	for block in d.block_start:
		intervals.append((block, block+2000))

	if len(intervals) > 0:
		putative_spots = list(find_intervals(intervals))
		chr = d.chr.tolist()[0]
		
		for (x, y) in putative_spots:
			if y - x < 5000:
				flank_rate = d[(d.block_start >= x) & (d.block_end <= y)].flank_rate.mean()
				spot_heat = d[(d.block_start >= x) & (d.block_end <= y)].rate_ratio.mean()
				
				o.write('%s,%s,%s,%.4f,%.4f\n'  % (chr, x, y, flank_rate, spot_heat))
o.close()
