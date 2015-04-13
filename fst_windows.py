import gzip
import re
import argparse
import numpy as np
import os

parser = argparse.ArgumentParser()
parser.add_argument("--comp", help="comparison for which to run analysis")
parser.add_argument("--chr", help="chromosome for which to run analysis")
args = parser.parse_args()
comp = args.comp
chr = args.chr

window = 10000

chr_lengths = { 'chr10': 20806668, 'chr11': 21403021, 'chr12': 21576510, 'chr13': 16962381,
                'chr14': 16419078, 'chr15': 14428146, 'chr16': 9909, 'chr17': 11648728,
                'chr18': 11201131, 'chr19': 11587733, 'chr1A': 73657157, 'chr1B': 1083483,
                'chr1': 118548696, 'chr20': 15652063, 'chr21': 5979137, 'chr22': 3370227,
                'chr23': 6196912, 'chr24': 8021379, 'chr25': 1275379, 'chr26': 4907541,
                'chr27': 4618897, 'chr28': 4963201, 'chr2': 156412533, 'chr3': 112617285,
                'chr4A': 20704505, 'chr4': 69780378, 'chr5': 62374962, 'chr6': 36305782,
                'chr7': 39844632, 'chr8': 27993427, 'chr9': 27241186, 'chrLG2': 109741,
                'chrLG5': 16416, 'chrLGE22': 883365, 'chrZ': 72861351}


def get_fst(comp, chr, window):
        fst_file = '/mnt/gluster/home/sonal.singhal1/LTF/analysis/fst/%s.Wrights_Fst.%s.csv.gz' % (chr, comp)
        snps = {}
	if os.path.isfile(fst_file):
                f = gzip.open(fst_file, 'r')
                header = f.next()
                for l in f:
                        d = re.split(',', l.rstrip())
                        fst = d[2]
                        if fst != 'nan':
                                if fst < 0:
                                        fst = 0
				snp_window = int(int(d[1])/float(window))
				if snp_window not in snps:
					snps[snp_window] = []
                                snps[snp_window].append(float(fst))
                f.close()
	return snps

def calc_windows(out, chr, chr_length, window, snps):
	o = open(out, 'w')
	o.write('chr,start,end,num_snps,avg_fst\n')
	for ix, start in enumerate(range(1, chr_length, window)):
		end = start + window - 1
		if end > chr_length:
			end = chr_length

		num_snps = 0
		avg_fst = 'nan'
		if ix in snps:
			num_snps = len(snps[ix])
			avg_fst = np.mean(snps[ix])

		o.write('%s,%s,%s,%s,%s\n' % (chr, start, end, num_snps, avg_fst))
	o.close()

snps = get_fst(comp, chr, window)
out = '/mnt/gluster/home/sonal.singhal1/LTF/analysis/fst/%s.Wrights_Fst.%s.window%s.csv' % (chr, comp, window)
calc_windows(out, chr, chr_lengths[chr], window, snps)

			


	
