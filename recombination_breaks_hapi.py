import re
from itertools import izip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--chr", help="chromosome for which to run analysis")
args = parser.parse_args()
chr = args.chr

dir = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/'
out = '%sbreaks/recombinationbreaks.%s.hapi.noswitch.csv' % (dir, chr)
ditch = '%sbreaks/badintervals.%s.hapi.noswitch.csv' % (dir, chr)

o = open(out, 'w')
o.write('chromosome,breaks,big_breaks,all_sites,het_sites\n')

sites_file = '%sinput_files/%s.hapi.noswitch.locs' % (dir, chr)
hap_file = '%soutput_files/%s.noswitch.csv' % (dir, chr)

# get sites
sites = []
f = open(sites_file, 'r')
for l in f:
	d = re.split('\s', l)
	sites.append(int(d[2]))
f.close()

# store haplotypes
haps = {}
f = open(hap_file, 'r')
for l in f:
	d = re.split(',', l.rstrip())
	hap = d[1:-1]
	if d[0] not in haps:
		haps[d[0]] = []
	haps[d[0]].append(hap)
f.close()

children = ['21', '22', '23']
parents = ['12', '10']

def delete_bp(bp, val1, val2, matches, locs, bad_breaks):
	p1 = re.compile('%s%s{1,%s}%s' % (val1, val2, bp, val1))
	spans = [x.span() for x in p1.finditer(matches)]
	# need to do this because deleting a list by index
	# so always need to start with backend to front
	spans.reverse()
	for span in spans:
		s = list(span)
		s[1] = s[1] - 1
		bad_breaks.append((locs[s[0]], locs[s[1]]))
		# again, deleting indices, so need to go from back to front
		for i in reversed(range(s[0] + 1, s[1])):
			del locs[i]
	matches = re.sub('%s%s{1,%s}%s' % (val1, val2, bp, val1), '%s%s' % (val1, val1), matches)
	return matches, locs, bad_breaks

# borrowed from 
# http://nbviewer.ipython.org/github/goldenhelix/cftr-variant-classification-analysis/blob/master/CFTR%20Data%20Munging.ipynb
def find_intervals(intervals):
	cur_start, cur_stop = intervals[0]
	for next_start, next_stop in intervals[1:]:
		if cur_stop <= next_start:
			yield cur_start, cur_stop
			if next_stop > cur_stop:
				cur_start, cur_stop = next_start, next_stop
			else:
				cur_start = next_start
        	else:
			if next_stop > cur_stop:
            			cur_stop = next_stop
	yield cur_start, cur_stop

bad_breaks = []
# identify heterozygous sites in child
for child in children:
	# paternal haplotype always comes first
	for parent, child_hap_wmiss in zip(parents, haps[child]):	
		match = []
		locs = []
		for ix, (a, b, c) in enumerate(izip(child_hap_wmiss, haps[parent][0], haps[parent][1])):
			if a != '0' and b != '0' and c != '0':
				if b != c:
					locs.append(sites[ix])
					if a == b:
						match.append('0')
					elif a == c:
						match.append('1')
		
		matches = ''.join(match)
		for bp in range(1,11):
			matches, locs, bad_breaks = delete_bp(bp, '0', '1', matches, locs, bad_breaks)
			matches, locs, bad_breaks = delete_bp(bp, '1', '0', matches, locs, bad_breaks)

		switch_all = 0
		switch_all_locs = []
		for ix, (pos1, pos2) in enumerate(zip(match, match[1:])):
			if pos1 != pos2:
				switch_all += 1
				
		switch_big = 0
		switch_big_locs = []
		for ix, (pos1, pos2) in enumerate(zip(matches, matches[1:])):
			if pos1 != pos2:
				switch_big += 1
		
		o.write('%s,%s,%s,%s,%s\n' % (chr, switch_all, switch_big, len(child_hap_wmiss), len(match)))
o.close()

d = open(ditch, 'w')
d.write('interval_start,interval_end\n')
if len(bad_breaks) > 1:
	breaks = list(find_intervals(sorted(bad_breaks)))
else:
	breaks = bad_breaks
for a, b in breaks:
	d.write('%s,%s\n' % (a, b))
d.close()
