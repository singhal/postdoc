import pandas as pd
import re
import numpy as np
from itertools import izip

gff = '/mnt/gluster/home/sonal.singhal1/reference/Taeniopygia_guttata.gff'
# long autosomal chromosomes
chrs = ['13', '11', '5', '7', '6', '1A', '4', '3', '2', '1', '8']
file_base = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/maps/chr%s_recombination_bpen100.rm.txt'
out = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/TSS/rho_near_TSS_1000kb.csv'
# plot 10^scale Kb
scale = 3

chr_lengths = { '10': 20806668, '11': 21403021, '12': 21576510, '13': 16962381,
                '14': 16419078, '15': 14428146, '16': 9909, '17': 11648728,
                '18': 11201131, '19': 11587733, '1A': 73657157, '1B': 1083483,
                '1': 118548696, '20': 15652063, '21': 5979137, '22': 3370227,
                '23': 6196912, '24': 8021379, '25': 1275379, '26': 4907541,
                '27': 4618897, '28': 4963201, '2': 156412533, '3': 112617285,
                '4A': 20704505, '4': 69780378, '5': 62374962, '6': 36305782,
                '7': 39844632, '8': 27993427, '9': 27241186, 'LG2': 109741,
                'LG5': 16416, 'LGE22': 883365, 'Z': 72861351}


###############

d = pd.read_csv(gff, sep='\t', header=0, names=['chr', 'type', 'cds_mrna', 
						'start', 'stop', 'score', 
						'orientation', 'codon_pos', 'id'])
gdata = d[ d.cds_mrna == 'mRNA' ].set_index('id').to_dict()
genes = d[ d.cds_mrna == 'mRNA' ].id.unique()

tss = {}
for gene in genes:
	chr = gdata['chr'][gene]
	orientation = 1
	if chr not in tss:
		tss[ chr ] = {}
	if gdata['orientation'][gene] == '+':
		tss_start = gdata['start'][gene]
	elif gdata['orientation'][gene]  == '-':
		tss_start = gdata['stop'][gene]
		orientation = -1
	tss[ chr ][ tss_start ] = orientation
	
negative = [-1 * x for x in np.logspace(0, scale, num=30)]
negative.reverse()
bins =  negative + [0] + np.logspace(0, scale, num=30).tolist()
bins = [ int(x * 1000) for x in bins]

o = open(out, 'w')
o.write('chr,tss,leftbin,rightbin,mean_tss_rho,mean_bin_rho\n')
for chr in chrs:
	if chr in tss:
		file = file_base % chr
		d = pd.read_csv(file, sep=" ", skiprows=3, header=None, 
			names=['left_snp', 'right_snp', 'meanrho', 'p0.05', 'p0.95'])
		for start in tss[chr]:
			tmp_bins = [start + x for x in bins]
			tmp_d = d[ d.right_snp >= tmp_bins[0] ]
			tmp_d = tmp_d[ tmp_d.right_snp <= tmp_bins[-1] ]

			global_rho = sum((tmp_d.right_snp - tmp_d.left_snp) * tmp_d.meanrho) / \
					sum((tmp_d.right_snp - tmp_d.left_snp))
			
			groups = tmp_d.groupby(np.digitize(tmp_d.right_snp, tmp_bins))

			for ix, group in groups:
				mean_rho =  sum((group.right_snp - group.left_snp) * group.meanrho) / \
						sum((group.right_snp - group.left_snp))

				if ix < len(bins):
					leftbin = bins[ix - 1] * tss[chr][start]
					rightbin = bins[ix] * tss[chr][start]
					
					rangebin = sorted([leftbin, rightbin])
					
					o.write('%s,%s,%s,%s,%s,%s\n' % 
							(chr, start, rangebin[0], rangebin[1], global_rho, mean_rho)) 				
o.close()

