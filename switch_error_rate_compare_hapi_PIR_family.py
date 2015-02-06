import re
import argparse
import gzip
from itertools import izip
import numpy as np

chr_lengths = { 'chr10': 20806668, 'chr11': 21403021, 'chr12': 21576510, 'chr13': 16962381,
                'chr14': 16419078, 'chr15': 14428146, 'chr16': 9909, 'chr17': 11648728,
                'chr18': 11201131, 'chr19': 11587733, 'chr1A': 73657157, 'chr1B': 1083483,
                'chr1': 118548696, 'chr20': 15652063, 'chr21': 5979137, 'chr22': 3370227,
                'chr23': 6196912, 'chr24': 8021379, 'chr25': 1275379, 'chr26': 4907541,
                'chr27': 4618897, 'chr28': 4963201, 'chr2': 156412533, 'chr3': 112617285,
                'chr4A': 20704505, 'chr4': 69780378, 'chr5': 62374962, 'chr6': 36305782,
                'chr7': 39844632, 'chr8': 27993427, 'chr9': 27241186, 'chrLG2': 109741,
                'chrLG5': 16416, 'chrLGE22': 883365}

parser = argparse.ArgumentParser(description='Identify switch errors between family samples phased via HAPI and phased via ShapeIT PIR.')
parser.add_argument('--chr', help='chromosome')
args = parser.parse_args()
chr = args.chr

dir = '/mnt/gluster/home/sonal.singhal1/ZF/'
# PIR output
pir_out = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/finch21/%s_haplotypes.haps' % chr
# HAPI output
hapi_out = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/output_files/%s.noswitch.csv' % chr
loc_out = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/input_files/%s.hapi.noswitch.locs' % chr
# out
out = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/switch_error/%s_switches.csv' % chr

# identify the sites phased via hapi
hapi_sites = {}
f = open(loc_out, 'r')
for l in f:
	d = re.split('\s+', l.rstrip())
	hapi_sites[int(d[2])] = 1
f.close()

# get sites phased by PIR
pir_sites = {}
pir_haps = {}
p = open(pir_out, 'r')
for l in p:
	d = re.split('\s+', l.rstrip())
	if int(d[2]) in hapi_sites:
		# this is a singleton :(
		# don't want to use singletons ...
		if d[5:].count('0') == 1 or d[5:].count('1') == 1:
			pass
		else:
			pir_sites[int(d[2])] = 1
			# the individuals are inverse in PIR and hapi, so need to invert them here
			for ix, site in zip([3,2,1,0], d[-4:]):
				if ix not in pir_haps:
					pir_haps[ix] = list()
				pir_haps[ix].append(site)
p.close()

hapi_sites_sorted = sorted(hapi_sites.keys())
hapi_haps = {}
h = open(hapi_out, 'r')
sites = []
for ix, l in enumerate(h):
	if ix < 4:
		d = re.split(',', l.rstrip())[1:-1]
		if ix not in hapi_haps:
			hapi_haps[ix] = []
		for site, pos in zip(hapi_sites_sorted, d):
			if site in pir_sites:
				sites.append(site)
				hapi_haps[ix].append(str(int(pos) - 1))
	else:
		pass
h.close()

def calculate_switch(ix1, ix2, ix3):
	matches = []
	het_sites = []
	for site, a1, a2, b in izip(sites, hapi_haps[ix1], hapi_haps[ix2], pir_haps[ix3]):
		if a1 == '-1' or a2 == '-1' or b == '-1':
			matches.append('-')
			het_sites.append(site)
		else:
			if a1 != a2:
				het_sites.append(site)
				if a1 != b:
					matches.append('F')
				else:
					matches.append('T')
	
	switch_locs = []
	for site, pos1, pos2 in izip(het_sites, matches, matches[1:]):
		if pos1 != '-' and pos2 != '-':
			if pos1 != pos2:
				switch_locs.append(site)
		
	return switch_locs, het_sites

switch_locs1, het_sites1 = calculate_switch(0,1,0)
switch_locs2, het_sites2 = calculate_switch(2,3,2)

bins = range(0,chr_lengths[chr],50000)
if len(bins) == 1:
	bins.append(50000)
site_bins = np.digitize(sites, bins)
het_bins1 = np.digitize(het_sites1, bins)
het_bins2 = np.digitize(het_sites2, bins)
switch_bins1 = []
switch_bins2 = []
if len(switch_locs1) > 0:
	switch_bins1 = np.digitize(switch_locs1, bins)
if len(switch_locs2) > 0:
	switch_bins2 = np.digitize(switch_locs2, bins)

o = open(out, 'w')
o.write('chromosome,start,end,number_snps,number_het_sites1,number_switches1,number_het_sites2,number_switches2\n')
for i in range(1, len(bins) + 1):
	bin_start = 1 + 50000 * (i - 1)
	bin_end = 50000 * i
	if i + 1 == len(bins) + 1:
		bin_end = chr_lengths[chr]
	
	o.write('%s,%s,%s,%s,%s,%s,%s,%s\n' % (chr, bin_start, bin_end, list(site_bins).count(i), \
						list(het_bins1).count(i), list(switch_bins1).count(i), \
						list(het_bins2).count(i), list(switch_bins2).count(i)))


