import re
import glob
import numpy as np
import pandas as pd
import random
from itertools import izip
import gzip
import os

putative_hotspots = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/seqldhot_hotspots/spot2kb_flank40kb.seqldhot_validate_hotspots.csv'
chrs = ['chr1', 'chr2', 'chr3']

vcfbase = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.%s.coverage.repeatmasked.filtered.nomendel.phased.vcf.gz'
rho_dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/without_fam/maps/'

# take this much sequence around the putative hotspot
block = 50e3
results_dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/hotspot_checks/seqldhot_sensitivity/'
theta = 0.00675

# get the most conservative set of hotspots 
d = pd.read_csv(putative_hotspots)
d = d[d.species == 'ZF']
d = d[d.zmatch_lk >= 10]
d = d[d.chr.isin(chrs)]

sh = open('%scommands.sh' % results_dir, 'w')
# now start to create the files
for chr, chrhot in d.groupby('chr'):
	vcf_file = vcfbase % chr

	rho_file = '%s%s.window10000.bpen100.rm.txt' % (rho_dir, chr)
	rhos = pd.read_csv(rho_file)
	haplo = {}

	f = gzip.open(vcf_file, 'r')
	for l in f:
		if not re.search('#', l):
			d = re.split('\s+', l.rstrip())

			# phased site
			if re.search('\|', ' '.join(d[9:])):
				pos = int(d[1])

				keep = True
				genos = []

				# multiallelic site
				if re.search(',', d[4]):
					keep = False

				for geno in d[9:]:
					geno = re.search('^([^:]+)', geno).group(1)
					genos = genos + re.split('\|', geno)

				count0 = genos.count('0')
				count1 = genos.count('1')

				# singleton
				if count0 ==  1 or count1 == 1:
					keep = False

				# these are the only sites we want
				# they are within the chunk to analyze, are greater than singletons and aren't in rm
				if keep:
					for ix, base in enumerate(genos):
						if ix not in haplo:
							haplo[ix] = {}
						haplo[ix][pos] = base
	f.close()

	for start in chrhot.spot_start:
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
	
			for ix, prop in enumerate([0.5,1.5]):
				results = '%sputative_hotspot_%s_%s.%s.seqLDhot' % (results_dir, chr, start, prop)
				out_in = '%sputative_hotspot_%s_%s.%s.seqLDhot.in' % (results_dir, chr, start, prop)
				o = open(out_in, 'w')
				o.write('Number of runs = 5000\nMIN number of iterations per hotspots = 100\ndriving values (for rho) = 2\n')
				o.write('background rho = %.3f\ntheta (per site) = %s\n' % (back_rho * 1000 * prop, theta))
				o.write('abs grid for hotspot likelihood\n0.5 40\nrel grid for hotspots likelihood\n5 100\n')
				o.write(' sub-region (number of SNPS; length (bps); frequency (bps))\n')
				o.write(' 10 2000 1000\n#\n')
				o.close()
				sh.write('/mnt/lustre/home/sonal.singhal1/bin/sequenceLDhot/sequenceLDhot %s %s %s\n' % (out_in, out, results))
sh.close()
