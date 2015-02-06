import re
import pandas as pd
from itertools import izip
import gzip
import os
import numpy as np
import glob
import random

phase = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/phase_hotspots/spot2kb_flank40kb.phase_validate_hotspots.csv'
seqld = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/seqldhot_hotspots/spot2kb_flank40kb.seqldhot_validate_hotspots.csv'
results_dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/hotspot_checks/phase_sensitivity/'
vcfbase = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.%s.coverage.repeatmasked.filtered.nomendel.phased.vcf.gz'
block = 50e3
theta = 0.00675
rho_dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/without_fam/maps/'

chrs = ['chr1', 'chr1A', 'chr2', 'chr3', 'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', \
                'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15']

valid = {}
# get the hotspots validated by seqld
d = pd.read_csv(seqld)
d = d[d.species == 'ZF']
d = d[d.zmatch_lk > 10]
for chr, start, rel_start, length in izip(d.chr, d.spot_start, d.zmatch_start, d.zmatch_length):
	rel_start = start + (rel_start - 25000)
	if chr not in valid:
		valid[chr] = dict()
	valid[chr][start] = {'rel_start': rel_start, 'length': length}

d = pd.read_csv(phase)
d = d[d.species == 'ZF']
d = d[d.zposterior < 0.9]
d = d[d.chr.isin(chrs)]

sh = open('%scommands.sh' % results_dir, 'w')
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

	for start in chrhot.spot_start:
		if start in valid[chr]:
			real_start = valid[chr][start]['rel_start']
			out = '%sputative_hotspot_%s_%s.PHASE.txt' % (results_dir, chr, start)
			if not os.path.isfile(out):
				seq_start = real_start - ( block / 2.0 )
				if seq_start < 1:
					seq_start = 1
				seq_end = real_start + ( block / 2.0 )

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
                     
                        	start_hot = int(block / 2.0) - 50
                        	stop_hot = start_hot + valid[chr][start]['length'] + 50
                        	sh.write('/mnt/lustre/home/sonal.singhal1/bin/phase.2.1.1/PHASE -X100 -k999 -MR2 1 %s %s -r%s %s %s 100 1 100\n' % \
					(start_hot, stop_hot, out_in, out, out2))
sh.close()

									
