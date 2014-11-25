import re
import pandas as pd
import os
from itertools import izip

putative_hotspots = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/LTF_ZF.putative_hotspots.csv'
chrs = ['chr1', 'chr1A', 'chr2', 'chr3', 'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', \
                'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chrZ']
out = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/seqldhot_hotspots/spot2kb_flank40kb.seqldhot_validate_hotspots.csv'

lr_cutoff = 10
center = 25000
dist = 2000


# get the most conservative set of hotspots 
d = pd.read_csv(putative_hotspots)
d = d[(d.spot_size == 2000) & (d.flank_size == 40000)]
d = d[d.chr.isin(chrs)]
d['rho_lambda'] = d.block_rate / d.flank_rate

hotspots = d.to_dict()
indices = hotspots['chr'].keys()

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

def parse_file(file, lr_cutoff, center):
	match_heat = None
	match_start = None
	match_length = None
	match_lk = None

	d = pd.read_csv(file, sep='\s+', skiprows=1, header=None)
	back_rho = d.X3.min()

	d_max = d[d.X2 >= lr_cutoff]
	intervals = [(i,j) for i,j in zip(d_max.X0, d_max.X1)]

	spots = []
	spot = []
	if len(intervals) > 0:
		if len(intervals) > 1:
			putative_spots = list(find_intervals(intervals))
		else:
			putative_spots = intervals
		for i, j in putative_spots:
			if abs(center - i) < dist:
				spots.append((i,j))
			else:
				if i < center and j > center:
					spots.append((i,j))
		if len(spots) > 1:
			heat = 0
			for tmpspot in spots:
				tmp = d[d.X0 >= tmpspot[0]]
				tmp = tmp[tmp.X1 <= tmpspot[1]]
				if tmp.X3.mean() > heat:
					heat = tmp.X3.mean()
					spot = tmpspot
		elif len(spots) == 1:
			spot = spots[0]

	if spot:
		d = d[d.X0 >= spot[0]]
		d = d[d.X1 <= spot[1]]
		match_lk = d.X2.mean()
		match_heat = d.X3.mean() / back_rho
		match_start = spot[0]
		match_length = spot[1] - spot[0]
					
	return (match_lk, match_heat, match_start, match_length)

o = open(out, 'w')
o.write(','.join(sorted(hotspots.keys())) + 'zmatch_lk,zmatch_heat,zmatch_start,zmatch_length,lmatch_lk,lmatch_heat,lmatch_start,lmatch_length\n')	
for index in indices:
	start = hotspots['spot_start'][index]
	chr = hotspots['chr'][index]
	rho_lambda = hotspots['rho_lambda'][index]
	block_rate = hotspots['block_rate'][index]

	out_zf = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/seqldhot_hotspots/putative_hotspot_%s_%s.seqLDhot.txt.sum' % (chr, start)
	out_ltf = '/mnt/gluster/home/sonal.singhal1/LTF/analysis/seqldhot_hotspots/putative_hotspot_%s_%s.seqLDhot.txt.sum' % (chr, start)	

	zmatch_lk, zmatch_heat, zmatch_start, zmatch_length = parse_file(out_zf, lr_cutoff, center)
	lmatch_lk, lmatch_heat, lmatch_start, lmatch_length = parse_file(out_ltf, lr_cutoff, center)
	o.write(','.join([str(hotspots[key][index]) for key in sorted(hotspots)]))
	o.write(',%s,%s,%s,%s,%s,%s,%s,%s\n' % (zmatch_lk, zmatch_heat, zmatch_start, \
			zmatch_length, lmatch_lk, lmatch_heat, lmatch_start, lmatch_length))
o.close()
		
		
