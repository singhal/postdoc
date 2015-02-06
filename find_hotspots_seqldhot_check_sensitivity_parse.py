import re
import glob
import pandas as pd
import os
from itertools import izip

files = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/hotspot_checks/seqldhot_sensitivity/*sum')

lr_cutoff = 10
center = 25000
dist = 5000

# find intervals of sequences that satisfy cut-offs
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

print 'chr,spot_start,back_rho_rate,match_lk,match_heat,match_start,match_length'
for file in files:
	chr = re.search('(chr[0-9|A-Z]+)', file).group(1)
	bp = re.search('chr[0-9|A-Z]+_(\d+)', file).group(1)
	rho = re.search('(\d\.\d)\.seqLDhot', file).group(1)

	match_lk, match_heat, match_start, match_length = parse_file(file, lr_cutoff, center)
	if isinstance(match_start, int):
		match_start = int(bp) + (match_start - center)
	print '%s,%s,%s,%s,%s,%s,%s' % (chr, bp, rho, match_lk, match_heat, match_start, match_length)
	
