import glob
import re
import pandas as pd
import os
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--sp", help="species for which to run the analysis")
args = parser.parse_args()
sp = args.sp

# repeat masked file because the samples aren't repeat masked
repeat = '/mnt/gluster/home/sonal.singhal1/reference/taeGut1.repeatMaskerBlast.repeatLibrary20140131.out'
# to get the background rate
rho_dir = '/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/maps/' % sp
# take this much sequence around the putative hotspot
block = 50e3
# output directory
results_dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/hotspot_checks/phasing_uncertainty/'
# command sh
command_sh = results_dir + 'seqldhot.sh'
# number of samples to test
samp = 3
# file with validated hotspots
d = pd.read_csv('/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/seqldhot_hotspots/spot2kb_flank40kb.seqldhot_validate_hotspots.csv')
spots = {}

chrs = [    'chr1', 'chr1A', 'chr2', 'chr3', 'chr4', 'chr4A', \
            'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', \
            'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chrZ']

if sp == 'LTF':
        theta = 0.0046
	d = d[ d.lmatch_lk >= 10]
	for chr, start in zip(d.chr, d.lmatch_start):
                if chr not in spots:
                        spots[chr] = []
                spots[chr].append(start)
if sp == 'ZF': 
        theta = 0.0058
	d = d[ d.zmatch_lk >= 10]
        for chr, start in zip(d.chr, d.zmatch_start):
		if chr not in spots:
			spots[chr] = []
		spots[chr].append(start)

# store all repeats into a dictionary
repeats = {}
f = open(repeat, 'r')
for l in f:
        if re.search('^\s+\d+', l):
                d = re.split('\s+', l.rstrip())
                if d[5] in chrs:
                        if d[5] not in repeats:
                                repeats[d[5]] = {}
                        start = int(d[6])
                        end = int(d[7]) + 1

                        for pos in range(start, end):
                                repeats[d[5]][pos] = 1
f.close()

sh = open(command_sh, 'w')

for chr in spots:
	rho_file = '%s%s.window10000.bpen100.rm.txt' % (rho_dir, chr)
	rhos = pd.read_csv(rho_file)
	
	for i in range(1,samp+1):
		haps = '/mnt/gluster/home/sonal.singhal1/%s/phasing/phasing_uncertainty/samples/%s.sample%s.haps' % (sp, chr, i)
		
		haplo = {}

		f = open(haps, 'r')
		for l in f:
			d = re.split('\s+', l.rstrip())
			# not a repeat
			if int(d[2]) not in repeats[chr]:
				count0 = d[5:].count('0')
				count1 = d[5:].count('1')
				# not a singleton
				if count0 > 1 and count1 > 1:
					for ix, base in enumerate(d[5:]):
						if ix not in haplo:
							haplo[ix] = {}
						haplo[ix][int(d[2])] = base
		f.close()

		for start in spots[chr]:
                	out = '%sputative_hotspot_%s_%s_%s.%s.seqLDhot.txt' % (results_dir, chr, start, i, sp)
			out_in = out.replace('txt', 'in')

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
                        	back_rho = round(np.mean(back_rho) * 1000, 3)
				# cannot handle back_rho of zero
				if back_rho == 0:
					back_rho = 0.001			

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

				o = open(out_in, 'w')
                        	o.write('Number of runs = 5000\nMIN number of iterations per hotspots = 100\ndriving values (for rho) = 2\n')
				o.write('background rho = %.3f\ntheta (per site) = %s\n' % (back_rho, theta))
                        	o.write('abs grid for hotspot likelihood\n0.5 40\nrel grid for hotspots likelihood\n5 100\n')
                        	o.write(' sub-region (number of SNPS; length (bps); frequency (bps))\n')
                        	o.write(' 10 2000 1000\n#\n')
                        	o.close()
                        	sh.write('/mnt/lustre/home/sonal.singhal1/bin/sequenceLDhot/sequenceLDhot %s %s\n' % (out_in, out))
sh.close()
