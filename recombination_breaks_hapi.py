import re
from itertools import izip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--chr", help="chromosome for which to run analysis")
args = parser.parse_args()
chr = args.chr

dir = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/'
out = '%srecombinationbreaks.%s.hapi.csv' % (dir, chr)

o = open(out, 'w')
o.write('chromosome,breaks,big_breaks,all_sites,het_sites,break_locations,big_break_locations\n')

sites_file = '%s%s.hapi.sites' % (dir, chr)
hap_file = '%shapi.%s.csv' % (dir, chr)

# identify the sites
sites = []
f = open(sites_file, 'r')
for l in f:
	d = re.split('\s+', l.rstrip())
	sites.append(float(d[2]))
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
		for bp in range(1,6):
			matches = re.sub('01{%s}0' % bp, '00', matches)
			matches = re.sub('10{%s}1' % bp, '11', matches)
				
		switch_all = 0
		switch_all_locs = []
		for ix, (pos1, pos2) in enumerate(zip(match, match[1:])):
			if pos1 != pos2:
				switch_all += 1
				switch_all_locs.append(locs[ix])
				
		switch_big = 0
		switch_big_locs = []
		for ix, (pos1, pos2) in enumerate(zip(matches, matches[1:])):
			if pos1 != pos2:
				switch_big += 1
				switch_big_locs.append(locs[ix])
		
		o.write('%s,%s,%s,%s,%s,%s,%s\n' % (chr, switch_all, switch_big, len(child_hap_wmiss), len(match), ':'.join([str(x) for x in switch_all_locs]), ':'.join([str(x) for x in switch_big_locs])))
o.close()
