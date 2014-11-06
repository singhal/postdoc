import re
import glob
import numpy as np
import pandas as pd
import random
from itertools import izip

putative_hotspots = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/LTF_ZF.putative_hotspots.csv'
chrs = ['chr1', 'chr1A', 'chr2', 'chr3', 'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', \
		'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15']
num_hotspots = 20
repeat_file = '/mnt/gluster/home/sonal.singhal1/reference/taeGut1.repeatMaskerBlast.repeatLibrary20140131.out'
hap_base = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/results/%s_haplotypes.haps'
# take this much sequence around the putative hotspot
block = 50e3
nhap = 38
results_dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/'
theta = 0.00654

# get a random set of hotspots 
d = pd.read_csv(putative_hotspots)
d = d[(d.species == 'ZF') & (d.spot_size == 2000) & (d.flank_size == 40000)]
d = d[d.chr.isin(chrs)]
d = d.ix[random.sample(d.index, num_hotspots)]

# remove repeat masked
repeats = {}
f = open(repeat_file, 'r')
for l in f:
        if re.search('^\s+\d+', l):
                a = re.split('\s+', l.rstrip())
                if a[5] in chrs:
                        if a[5] not in repeats:
                                repeats[a[5]] = {}
                        start = int(a[6])
                        end = int(a[7]) + 1

                        for pos in range(start, end):
                                repeats[a[5]][str(pos)] = 1
f.close()

for chr in chrs:
	repeats[chr] = {}

sh = open('%scommands.sh' % results_dir, 'w')
# now start to create the files
for chr, start, back_rho in zip(d.chr, d.spot_start, d.flank_rate):
	hap_file = hap_base % chr
	seq_start = start - ( block / 2.0 )
	if seq_start < 1:
		seq_start = 1
	seq_end = start + ( block / 2.0 )

	sites = {}
	haplo = {}

	for i in range(nhap):
		haplo[i] = []

	# in order to be able to make the input file for all, need:
	#	(1) to track the sites
	# 	(2) to track 0 or 1
	#	(3) to track ref and alt
	f = open(hap_file, 'r')
	for l in f:
		d = re.split('\s+', l.rstrip())
		
		pos = int(d[2])
		keep = True
		if pos < seq_start:
			keep = False
		if pos > seq_end:
			keep = False
		# not in repeat masked
		if pos in repeats[chr]:
			keep = False
		
		if keep:
			count0 = d[5:].count('0')
			count1 = d[5:].count('1')

			# singleton
			if count0 ==  1 or count1 == 1:
				keep = False

		# these are the only sites we want
		# they are within the chunk to analyze, are greater than singletons and aren't in rm
		if keep:
			sites[pos] = {'ref': d[3], 'alt': d[4]}
			for ix, base in enumerate(d[5:]):
				haplo[ix].append(int(base))
	f.close()

	# first, let's make the InferRho
	out = '%sputative_hotspot_%s_%s.InferRho.txt' % (results_dir, chr, start)
	o = open(out, 'w')
	o.write('1\nh\n\n%s\n%s\n' % (nhap, len(sites.keys())))
	sorted_sites = sorted(sites.keys())
	for ix in range(nhap):
		hap = ''
		for pos, allele in zip(sorted_sites, haplo[ix]):
			if allele == 1:
				hap += '%s ' % (sites[pos]['alt'])
			else:
				hap += '%s ' % (sites[pos]['ref'])
		o.write(hap + '\n')
	o.write(' '.join([str(x - int(seq_start) + 1) for x in sorted_sites]) + '\n')
	o.close()
	sh.write('/mnt/lustre/home/sonal.singhal1/bin/InferRho-1.0/InferRho %s\n' % out)

	# next, let's make the PHASE
	out = '%sputative_hotspot_%s_%s.PHASE.txt' % (results_dir, chr, start)
        o = open(out, 'w')
        o.write('%s\n%s\n' % (int(nhap/2.0), len(sites.keys())))
        sorted_sites = sorted(sites.keys())
	o.write('P ' + ' '.join([str(x - int(seq_start) + 1) for x in sorted_sites]) + '\n')
        o.write('S' * len(sorted_sites) + '\n')
	for ix in range(0,nhap,2):
                hap1 = ''.join([str(x) for x in haplo[ix]])
		hap2 = ''.join([str(x) for x in haplo[ix+1]])
		o.write('#%s\n%s\n%s\n' % (ix, hap1, hap2))
        o.close()
	out2 = out.replace('.txt', '.out')
	
	out_in = '%sputative_hotspot_%s_%s.PHASE.prior' % (results_dir, chr, start)	
	o = open(out_in, 'w')
	o.write('%s\n10\n1000\n200\n200\n50000\n' % back_rho)
	o.close()	
	
	sh.write('/mnt/lustre/home/sonal.singhal1/bin/phase.2.1.1/PHASE -k999 -MR1 1 -r%s %s %s 100 1 100\n' %(out_in,out,out2))

	# next, let's make the sequenceLDhot
        out = '%sputative_hotspot_%s_%s.seqLDhot.txt' % (results_dir, chr, start)
        o = open(out, 'w')
        o.write('Distinct = %s\nGenes = %s\nLoci = %s\n' % (nhap, nhap, len(sites.keys())))
	o.write('I=1\nK = -2\nPositions of loci:\n')
        sorted_sites = sorted(sites.keys())
        o.write(' '.join([str(x - int(seq_start) + 1) for x in sorted_sites]) + '\n')
        o.write('Haplotypes\n')
        for ix in range(nhap):
                hap = ''.join([str(x) for x in haplo[ix]])
                o.write('\t%s\t1\n' % (hap))
	o.write('#')
        o.close()
	
	out_in = '%sputative_hotspot_%s_%s.seqLDhot.in' % (results_dir, chr, start)
	o = open(out_in, 'w')
	o.write('Number of runs = 5000\nMIN number of iterations per hotspots = 100\ndriving values (for rho) = 2\n')
	o.write('background rho = %s\ntheta (per site) = %s\n' % (back_rho * 1000, theta))
	o.write('abs grid for hotspot likelihood\n0.5 40\nrel grid for hotspots likelihood\n5 100\n')
	o.write(' sub-region (number of SNPS; length (bps); frequency (bps))\n')
	o.write(' 10 2000 1000\n#\n')
	o.close()
	sh.write('/mnt/lustre/home/sonal.singhal1/bin/sequenceLDhot/sequenceLDhot %s %s\n' % (out_in, out))
sh.close()
