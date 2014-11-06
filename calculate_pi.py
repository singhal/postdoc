import re
import gzip
import glob
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--sp", help="species [ZF|LTF]")
parser.add_argument("--nind", help="number of haplotypes")

args = parser.parse_args()
sp = args.sp
nind = int(args.nind)

chr_lengths = {	'chr10': 20806668, 'chr11': 21403021, 'chr12': 21576510, 'chr13': 16962381,
                'chr14': 16419078, 'chr15': 14428146, 'chr16': 9909, 'chr17': 11648728,
                'chr18': 11201131, 'chr19': 11587733, 'chr1A': 73657157, 'chr1B': 1083483,
                'chr1': 118548696, 'chr20': 15652063, 'chr21': 5979137, 'chr22': 3370227,
                'chr23': 6196912, 'chr24': 8021379, 'chr25': 1275379, 'chr26': 4907541,
                'chr27': 4618897, 'chr28': 4963201, 'chr2': 156412533, 'chr3': 112617285,
                'chr4A': 20704505, 'chr4': 69780378, 'chr5': 62374962, 'chr6': 36305782,
                'chr7': 39844632, 'chr8': 27993427, 'chr9': 27241186, 'chrLG2': 109741,
                'chrLG5': 16416, 'chrLGE22': 883365, 'chrZ': 72861351}
hap_dir = '/mnt/gluster/home/sonal.singhal1/gene_trees/chromosomes/%s' % sp
window = 50000

outfile = '/mnt/gluster/home/sonal.singhal1/%s/analysis/pop_gen/pi.csv' % sp
o = open(outfile, 'w')
o.write('chr,index,start,end,pi\n')

for chr, length in chr_lengths.items():
	seq_file = '%s_%s_haplotypes.fasta' % (hap_dir, chr)
	for ix, start in enumerate(range(1, length, window)):
		end = start + window
		if end > length:
			end = length
		
		haps = []
		for i in range(nind):
			haps.append('')
			seq = subprocess.Popen('samtools faidx %s haplo%s:%s-%s' % (seq_file, i, start, end), shell=True, stdout=subprocess.PIPE)
			for l in seq.stdout:
				if not re.match('>', l):
					haps[i] += l.rstrip()
		
		pi_sum = 0
		for i in range(nind):
			for j in range(i):
				diff = 0
				seq_length = 0
				for bp1, bp2 in zip(haps[i], haps[j]):
					if bp1 != 'N' or bp2 != 'N':
						seq_length += 1
						if bp1 != bp2:
							diff += 1
				if seq_length > 0:
					pi_sum += (1/float(nind)) * (1/float(nind)) * (diff/float(seq_length))
		if seq_length > 0:
			pi_sum = pi_sum * 2
			o.write('%s,%s,%s,%s,%.4f\n' % (chr, ix, start, end, pi_sum))
		else:
			o.write('%s,%s,%s,%s,NA\n' % (chr, ix, start, end))
o.close()		
