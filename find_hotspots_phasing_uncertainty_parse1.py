import re
import glob
import pandas as pd
import os
from itertools import izip

files = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/hotspot_checks/phasing_uncertainty/*txt')

lr_cutoff = 10
center = 25000
dist = 3000

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


def parse_file(txt, file, lr_cutoff, center):
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
                match_start = spot[0]
                match_length = spot[1] - spot[0]
                                        
        return (match_lk, match_heat, match_start, match_length)

def check_file(file, out):
	good = False

	f = open(file, 'r')
	for l in f:
		if re.search('Positions', l):
			positions = f.next().rstrip()
			sites = [int(x) for x in re.split('\s+', positions)]
	f.close()

	last_snp = 0
	f = open(out, 'r')
	for l in f:
		if re.search('^\d', l):
			d = re.split('\s+', l.rstrip())
			if len(d) > 4:
				last_snp = int(d[-1])
	f.close()

	if len(sites) > 0:
		if abs(sites[-1] - last_snp) < 1000:
			good = True

	return good

print 'chr,putative_start,rep,lk,heat,start,length'
for file in files:
	chr = re.search('(chr[0-9|A-Z]+)', file).group(1)
	bp = re.search('chr[0-9|A-Z]+_(\d+)', file).group(1)
	rep = re.search('\.\d_(\d)\.', file).group(1)
	
	out = file + '.sum'

	try: 
		if os.path.getsize(out) > 0:
			if check_file(file, out):
				match_lk, match_heat, match_start, match_length = parse_file(file, out, lr_cutoff, center)
				match_start = int(bp) + (match_start - center)
				print '%s,%s,%s,%s,%s,%s,%s' % (chr, bp, rep, match_lk, match_heat, match_start, match_length)
	except:
		pass

