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
                'chrLG5': 16416, 'chrLGE22': 883365, 'chrZ': 72861351}

parser = argparse.ArgumentParser(description='Identify switch errors between family samples phased via HAPI and phased via ShapeIT PIR.')
parser.add_argument('--chr', help='chromosome')
args = parser.parse_args()
chr = args.chr

dir = '/mnt/gluster/home/sonal.singhal1/ZF/'
vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/all_vcf/gatk.ug.all_zf.%s.coverage.filtered.repeatmasked.vqsr.vcf.gz' % chr
# PIR output
pir_out = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/phasing_uncertainty/switch_error/%s_haplotypes.haps' % chr
# HAPI output
hapi_out = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/output_files/hapi.%s.csv' % chr
# out
out = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/phasing_uncertainty/switch_error/%s_switches.csv' % chr

errorfile = dir + 'mendelian_errors/plink_results/all_zf.me.%s.mendel' % chr
sites_file = dir + 'mendelian_errors/all_zf.%s.map' % chr

sites_f = open(sites_file, 'r')
sites = {}
for l in sites_f:
	d = re.split('\s+', l)
	sites[d[1]] = int(d[3])
sites_f.close()
        
errors = {}
error_f = open(errorfile, 'r')
junk = error_f.next()
for l in error_f:
	d = re.split('\s+', l)
	errors[sites[d[4]]] = 1
error_f.close()
sites = {}

# identify the sites phased via hapi
hapi_sites = {}
v = gzip.open(vcf, 'r')
for l in v:
	if not re.search('^#', l):
		d = re.split('\t', l.rstrip())
		alleles = [d[3]] + re.split(',', d[4])
		indel = False
		for allele in alleles:
			if len(allele) > 1:
				indel = True
		if not indel and len(alleles) == 2:
			if int(d[1]) not in errors:
				all_missing = True
				for geno in d[-5:]:
					geno = re.search('^(\S/\S)', geno).group(1)
					if geno != './.':
						all_missing = False
				if not all_missing:
					hapi_sites[int(d[1])] = 1
v.close()

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
site_bins = np.digitize(sites, bins)
het_bins1 = np.digitize(het_sites1, bins)
het_bins2 = np.digitize(het_sites2, bins)
switch_bins1 = np.digitize(switch_locs1, bins)
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


