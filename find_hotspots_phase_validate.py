import re
import glob
import numpy as np
import pandas as pd
import random
from itertools import izip
import argparse
import gzip
import os

parser = argparse.ArgumentParser()
parser.add_argument("--sp", help="species for which to run the analysis")
args = parser.parse_args()
sp = args.sp

putative_hotspots = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/LTF_ZF.putative_hotspots.csv'
chrs = ['chr1', 'chr1A', 'chr2', 'chr3', 'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', \
		'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr27', 'chrZ']

if sp == 'ZF':
	vcfbase = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.%s.coverage.repeatmasked.filtered.nomendel.phased.vcf.gz'
if sp == 'LTF':
	vcfbase = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.%s.filtered.coverage.repeatmasked.vqsr.phased.vcf.gz'
rho_dir = '/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/maps/' % sp

# take this much sequence around the putative hotspot
block = 50e3
results_dir = '/mnt/gluster/home/sonal.singhal1/%s/analysis/phase_hotspots/' % sp

if sp == 'LTF':
	theta = 0.0046
if sp == 'ZF': 
	theta = 0.0058

# get the most conservative set of hotspots 
d = pd.read_csv(putative_hotspots)
d = d[(d.spot_size == 2000) & (d.flank_size == 40000)]
d = d[d.chr.isin(chrs)]

sh = open('%scommands.sh' % results_dir, 'w')
# now start to create the files
for chr, chrhot in d.groupby('chr'):
	vcf_file = vcfbase % chr
	if chr == 'chrZ' and sp == 'ZF':
		vcf_file = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.chrZ.coverage.repeatmasked.filtered.recodedsex.vqsr.phased.vcf.gz'
	if chr == 'chrZ' and sp == 'LTF':
		vcf_file = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.chrZ.filtered.coverage.repeatmasked.recodedsex.vqsr.phased.vcf.gz'
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
					genos.append(re.split('\|', geno))

				count0 = 0
				count1 = 0
				for geno in genos:
					count0 += geno.count('0')
					count1 += geno.count('1')

				# singleton
				if count0 ==  1 or count1 == 1:
					keep = False

				# these are the only sites we want
				if keep:
					for ix, geno in enumerate(genos):
						if ix not in haplo:
							haplo[ix] = {}
						haplo[ix][pos] = geno
	f.close()

	for start, length in zip(chrhot.spot_start, chrhot.length):
		out = '%sputative_hotspot_%s_%s.PHASE.txt' % (results_dir, chr, start)

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

			o = open(out, 'w')
        		o.write('%s\n%s\n' % (int(len(haplo)), len(sites)))
		        o.write('P ' + ' '.join([str(x - int(seq_start) + 1) for x in sorted_sites]) + '\n')
        		o.write('S' * len(sorted_sites) + '\n')
        		for ix in haplo:
				# diploid or haploid?
				diploid = True
				if len(haplo[ix][haplo[ix].keys()[0]]) == 1:
					diploid = False
				if diploid:
					hap1 = ''
					hap2 = ''
					for pos in sorted_sites:
						add1 = '?'
						add2 = '?'
						if pos in haplo[ix]:
							if haplo[ix][pos][0] != '.':
								add1 = haplo[ix][pos][0]
								add2 = haplo[ix][pos][1]	
        	        			hap1 += add1
						hap2 += add2
                			o.write('#%s\n%s\n%s\n' % (ix, hap1, hap2))
				else:
					hap1 = ''
					for pos in sorted_sites:
                                                add1 = '?'
                                                if pos in haplo[ix]:
                                                        if haplo[ix][pos][0] != '.':
								add1 = haplo[ix][pos][0]
						hap1 += add1
					o.write('#%s\n%s\n' % (ix, hap1))
        		o.close()
        		out2 = out.replace('.txt', '.out')

			out_in = '%sputative_hotspot_%s_%s.PHASE.prior' % (results_dir, chr, start)     
		        o = open(out_in, 'w')
        		o.write('%s\n10\n1000\n200\n200\n50000\n' % back_rho)
        		o.close() 

			out2 = out.replace('.txt', '.asphased.out')
		        out3 = out.replace('.txt', '.rephased.out')

			start_hot = int(block / 2.0) - 50
			stop_hot = start_hot + length + 50
			sh.write('/mnt/lustre/home/sonal.singhal1/bin/phase.2.1.1/PHASE -X100 -k999 -MR2 1 %s %s -r%s %s %s 100 1 100\n' % (start_hot, stop_hot, out_in, out, out2))
			sh.write('/mnt/lustre/home/sonal.singhal1/bin/phase.2.1.1/PHASE -X100 -MR2 1 %s %s -r%s %s %s 100 1 100\n' % (start_hot, stop_hot, out_in, out, out3))			

sh.close()
