import re
import pandas as pd
import os
from itertools import izip
import glob
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--match_dist", help="match_dist")
args = parser.parse_args()

files = glob.glob('/mnt/gluster/home/sonal.singhal1/simulations/shared/seqldhot/ZF/*sum')

lr_cutoff = 10
center = 25000
# spots have to be within 3000 of original location to be considered fit
fit_dist = 3000
# spots have to be within 3000 to be considered matching
match_dist = int(args.match_dist)
# spots within 3000 of each other will be dropped
max_dist = 3000

out =  '/mnt/gluster/home/sonal.singhal1/simulations/shared/seqldhot_validate_hotspots.match_dist%s.csv' % match_dist

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
			if abs(center - i) < fit_dist:
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
		match_mid = match_start + int(match_length / 2.0)
					
	return (match_lk, match_heat, match_start, match_mid, match_length, back_rho / 1000.)


def run_script(file, sp, chr, start, lr_cutoff, center, spots):
	try:
		if os.stat(file).st_size > 0:
                        match_lk, match_heat, match_start, match_mid, match_length, match_back_rho = parse_file(file, start, lr_cutoff, center)
	
			if chr not in spots[sp]:
		                spots[sp][chr] = {}

			if match_lk != 'NA':
        			spots[sp][chr][match_mid] = {'start': match_start, 'lk': match_lk, 'heat': match_heat, 'length': match_length, 'flank_rho': match_back_rho}
        except:
		pass

	return spots


chrs = {}
spots = {'ZF': {}, 'LTF': {}}
for out_zf in files:
	rho_zf = float(re.search('rho([0-9|\.]+)', out_zf).group(1)) 
	rho_ltf = rho_zf / 2.0
	out_ltf = out_zf.replace('ZF','LTF')
	out_ltf = out_ltf.replace(str(rho_zf), str(rho_ltf))
	
	chr = re.search('(rho.*diff\d+_\d+)', out_zf).group(1)
	start = re.search('_(\d+)\.seq', out_zf).group(1)
	start = int(start)
	if chr not in chrs:
		chrs[chr] = 1

	spots = run_script(out_zf, 'ZF', chr, start, lr_cutoff, center, spots)
	spots = run_script(out_ltf, 'LTF', chr, start, lr_cutoff, center, spots)


def removeCloseItems(items, itemDistance):
    if items:
        lastOutput = items[0]
        yield items[0]
        for currentItem in items[1:]:
            if ((currentItem - lastOutput) > itemDistance):
                lastOutput = currentItem
                yield currentItem

o = open(out, 'w')
o.write('chr,zmid,zlength,zbackrho,zheat,zlk,lmid,llength,lbackrho,lheat,llk\n')

for chr in chrs:
	# first get rid of spots too close to each other
	if chr in spots['ZF']:
		zf_mids = sorted(spots['ZF'][chr])
		zf_mids = [x for x in removeCloseItems(zf_mids, max_dist)]
	else:
		zf_mids = []
	if chr in spots['LTF']:
		ltf_mids = sorted(spots['LTF'][chr])
		ltf_mids = [x for x in removeCloseItems(ltf_mids, max_dist)]
	else:
		ltf_mids = []

	matches = {}
	rev_matches = {}
	for zf_mid in zf_mids:
		dists = [abs(x - zf_mid) for x in ltf_mids]
		if len(dists) > 1:
			min_dist = np.min(dists)
		elif len(dists) == 1:
			min_dist = dists[0]
		if len(dists) > 0:
			if min_dist < match_dist:
				matches[zf_mid] = ltf_mids[dists.index(min_dist)]
				rev_matches[ltf_mids[dists.index(min_dist)]] = zf_mid

	for zmid in zf_mids:
		if zmid in matches:
			lmid = matches[zmid]
			o.write('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (chr, zmid, spots['ZF'][chr][zmid]['length'], 
								spots['ZF'][chr][zmid]['flank_rho'], 
								spots['ZF'][chr][zmid]['heat'],
								spots['ZF'][chr][zmid]['lk'], lmid, spots['LTF'][chr][lmid]['length'], 
								spots['LTF'][chr][lmid]['flank_rho'],
								spots['LTF'][chr][lmid]['heat'], spots['LTF'][chr][lmid]['lk']))
		else:
			o.write('%s,%s,%s,%s,%s,%s,NA,NA,NA,NA,NA\n' % (chr, zmid, spots['ZF'][chr][zmid]['length'], 
								spots['ZF'][chr][zmid]['flank_rho'],
								spots['ZF'][chr][zmid]['heat'],
                                                                spots['ZF'][chr][zmid]['lk']))
	
	for lmid in ltf_mids:
		if lmid not in rev_matches:
			o.write('%s,NA,NA,NA,NA,NA,%s,%s,%s,%s,%s\n' % (chr, lmid, spots['LTF'][chr][lmid]['length'],
									spots['LTF'][chr][lmid]['flank_rho'],
                                                                	spots['LTF'][chr][lmid]['heat'], spots['LTF'][chr][lmid]['lk']))
