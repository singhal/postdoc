import re
import glob
import numpy as np
import pandas as pd
import random
from itertools import izip
import argparse
import gzip
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--sp", help="species for which to run the analysis")
args = parser.parse_args()
sp = args.sp

putative_hotspots = glob.glob('/mnt/gluster/home/sonal.singhal1/*/analysis/LDhelmet/*heat5*out')
chrs = ['chr1', 'chr1A', 'chr2', 'chr3', 'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', \
	'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chrZ']
results_dir = '/mnt/gluster/home/sonal.singhal1/%s/analysis/hotspots/seqldhot_hotspots/' % sp

if sp == 'ZF':
	hapbase = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/finch19/%s_haplotypes.haps'
	rho_dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/without_fam/maps/'
	theta = 0.00675
if sp == 'LTF':
	hapbase = '/mnt/gluster/home/sonal.singhal1/LTF/phasing/PIR_approach/%s_haplotypes.haps'
	rho_dir = '/mnt/gluster/home/sonal.singhal1/LTF/analysis/LDhelmet/old/maps/'
	theta = 0.00472

# take this much sequence around the putative hotspot
block = 50e3

# get the most conservative set of hotspots
hotspots = {}
for putative_hotspot in putative_hotspots: 
	d = pd.read_csv(putative_hotspot)
	d = d[(d.block == 2000) & (d.flank == 40000)]
	for chr, start in zip(d.chr, d.spot_start):
		if chr not in hotspots:
			hotspots[chr] = []
		if len(hotspots[chr]) > 5000:
			min_dist = np.min([abs(x - start) for x in hotspots[chr]])
			if min_dist > 0:
				hotspots[chr].append(start)
		else:
			hotspots[chr].append(start)


sh = open('/mnt/lustre/home/sonal.singhal1/scripts/seqldhot_commands_%s.sh' % sp, 'w')
# now start to create the files
for chr in hotspots:
	hap_file = hapbase % chr
	rho_file = '%s%s.window10000.bpen100.rm.txt' % (rho_dir, chr)

	rhos = pd.read_csv(rho_file)
	
	haplo = {}

	f = open(hap_file, 'r')
	for l in f:
		d = re.split('\s+', l.rstrip())
		genos = d[5:]
		count0 = genos.count('0')
		count1 = genos.count('1')

		# ditch singleton
		if count0 >  1 and count1 > 1:
			for ix, base in enumerate(genos):
				if ix not in haplo:
					haplo[ix] = {}
				haplo[ix][int(d[2])] = base
	f.close()

	for start in hotspots[chr]:
		out = '%sputative_hotspot_%s_%s.seqLDhot.txt' % (results_dir, chr, start)

		if not os.path.isfile(out):
			seq_start = start - ( block / 2.0 )
        		if seq_start < 1:
                		seq_start = 1
			seq_end = start + ( block / 2.0 )

			sites = filter(lambda x: x >= seq_start, haplo[0].keys())
			sites = filter(lambda x: x <= seq_end, sites)

			back_rho = rhos[ rhos.window_end >= seq_start ]
			back_rho = back_rho[ back_rho.window_start  <= seq_end ].rate
			back_rho = filter(lambda x: np.isfinite(x), back_rho)
			back_rho = np.mean(back_rho)
			if np.isfinite(back_rho):
				if back_rho == 0:
					back_rho = 0.0001
			else:
				back_rho = 0.0001

			sorted_sites = sorted(sites)
			# next, let's make the haplotypes
			tmp_haplo = {}
			for ix in sorted(haplo.keys()):
                                hap = ''
                                for pos in sorted_sites:
                                        hap += haplo[ix][pos]
				if hap not in tmp_haplo:
					tmp_haplo[hap] = 0
				tmp_haplo[hap] += 1

			# next, let's make the sequenceLDhot
        		o = open(out, 'w')
        		o.write('Distinct = %s\nGenes = %s\nLoci = %s\n' % (len(tmp_haplo), len(haplo), len(sites)))
			o.write('I=1\nK = -2\nPositions of loci:\n')
			o.write(' '.join([str(x - int(seq_start) + 1) for x in sorted_sites]) + '\n')
        		o.write('Haplotypes\n')
        		for hap, count in tmp_haplo.items():
				o.write('\t%s %s\n' % (hap, count))			
			o.write('#')
        		o.close()
	
			out_in = '%sputative_hotspot_%s_%s.seqLDhot.in' % (results_dir, chr, start)
			o = open(out_in, 'w')
			o.write('Number of runs = 5000\nMIN number of iterations per hotspots = 100\ndriving values (for rho) = 2\n')
			o.write('background rho = %.3f\ntheta (per site) = %s\n' % (back_rho * 1000, theta))
			o.write('abs grid for hotspot likelihood\n0.5 40\nrel grid for hotspots likelihood\n5 100\n')
			o.write(' sub-region (number of SNPS; length (bps); frequency (bps))\n')
			o.write(' 10 2000 1000\n#\n')
			o.close()
			sh.write('/mnt/lustre/home/sonal.singhal1/bin/sequenceLDhot/sequenceLDhot %s %s\n' % (out_in, out))
sh.close()
