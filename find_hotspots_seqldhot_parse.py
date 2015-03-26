import re
import pandas as pd
import os
from itertools import izip
import glob
import numpy as np

files = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/seqldhot_hotspots/*txt')
chrs = ['chr1', 'chr1A', 'chr2', 'chr3', 'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', \
        'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chrZ']
out = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/spot2kb_flank40kb.seqldhot_validate_hotspots.csv'

lr_cutoff = 10
center = 25000
# spots have to be within 3000 of original location to be considered fit
dist = 3000
# spots have to be within 3000 to be considered matching
# spots within 3000 of each other will be dropped
max_dist = 3000

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

def parse_file(file, start, lr_cutoff, center):
	match_heat = 'NA'
	match_start = 'NA'
	match_length = 'NA'
	match_lk = 'NA'

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
		match_start = start + (spot[0] - 25000)
		match_length = spot[1] - spot[0]
					
	return (match_lk, match_heat, match_start, match_length, back_rho / 1000.)


def run_script(file, sp, chr, start, lr_cutoff, center, spots):
	try:
		if os.stat(file).st_size > 0:
                        match_lk, match_heat, match_start, match_length, match_back_rho = parse_file(file, start, lr_cutoff, center)
	
			if chr not in spots[sp]:
		                spots[sp][chr] = {}

			if match_lk != 'NA':
        			spots[sp][chr][match_start] = {'lk': match_lk, 'heat': match_heat, 'length': match_length, 'flank_rho': match_back_rho}
        except:
		pass

	return spots


spots = {'ZF': {}, 'LTF': {}}
for out_zf in files:
	out_zf = out_zf + '.sum'
	out_ltf = out_zf.replace('ZF','LTF')
	
	chr = re.search('(chr[0-9|A-Z]+)', out_zf).group(1)
	start = re.search('_(\d+)\.seq', out_zf).group(1)
	start = int(start)

	spots = run_script(out_zf, 'ZF', chr, start, lr_cutoff, center, spots)
	spots = run_script(out_ltf, 'LTF', chr, start, lr_cutoff, center, spots)


o = open(out, 'w')
o.write('chr,zstart,zlength,zbackrho,zheat,zlk,lstart,llength,lbackrho,lheat,llk\n')

for chr in chrs:
	# zf first, and then ltf
	matches = {}

	zf_spots = sorted(spots['ZF'][chr].keys())
	ltf_spots = sorted(spots['LTF'][chr].keys())

	for spot in zf_spots:
		dists = [abs(spot - x) for x in ltf_spots]
		min_dist = np.min(dists)
		if min_dist <= max_dist:
			matches[spot] = ltf_spots[dists.index(min_dist)]
	
	# want to get rid of any hotspots that are too close to each other
	for a, b in izip(zf_spots, zf_spots[1:]):
		dist = b - a
		if dist < max_dist:
			if a not in matches and b not in matches:
				del spots['ZF'][chr][b]
			if a not in matches and b in matches:
	 			del spots['ZF'][chr][a]
			if a in matches and b in matches:
				del spots['ZF'][chr][b]
		
	rev_matches = {}
	for a, b in matches.items():
		rev_matches[b] = a

	for a, b in izip(ltf_spots, ltf_spots[1:]):
                dist = b - a
                if dist < max_dist:
                        if a not in rev_matches and b not in rev_matches:
				del spots['LTF'][chr][b]
                        if a not in rev_matches and b in rev_matches:
                                del spots['LTF'][chr][a]
                        if a in rev_matches and b in rev_matches:
                                del spots['LTF'][chr][b]
                                del matches[rev_matches[b]]
	
	rev_matches = {}
        for a, b in matches.items():
                rev_matches[b] = a

	for zstart in spots['ZF'][chr]:
		if zstart in matches:
			lstart = matches[zstart]
			o.write('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (chr, zstart, spots['ZF'][chr][zstart]['length'], 
								spots['ZF'][chr][zstart]['flank_rho'], 
								spots['ZF'][chr][zstart]['heat'],
								spots['ZF'][chr][zstart]['lk'], lstart, spots['LTF'][chr][lstart]['length'], 
								spots['LTF'][chr][lstart]['flank_rho'],
								spots['LTF'][chr][lstart]['heat'], spots['LTF'][chr][lstart]['lk']))
		else:
			o.write('%s,%s,%s,%s,%s,%s,NA,NA,NA,NA,NA\n' % (chr, zstart, spots['ZF'][chr][zstart]['length'], 
								spots['ZF'][chr][zstart]['flank_rho'],
								spots['ZF'][chr][zstart]['heat'],
                                                                spots['ZF'][chr][zstart]['lk']))
	
	for lstart in spots['LTF'][chr]:
		if lstart not in rev_matches:
			o.write('%s,NA,NA,NA,NA,NA,%s,%s,%s,%s,%s\n' % (chr, lstart, spots['LTF'][chr][lstart]['length'],
									spots['LTF'][chr][lstart]['flank_rho'],
                                                                spots['LTF'][chr][lstart]['heat'], spots['LTF'][chr][lstart]['lk']))
