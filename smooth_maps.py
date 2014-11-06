import pandas as pd
import argparse
import numpy as np
from itertools import izip

def find_hotspots(file, out, chr, block_size):
	d = pd.read_csv(file, sep=" ", skiprows=3, header=None, 
		names=['left_snp', 'right_snp', 'meanrho', 'p0.05', 'p0.95'])
	o = open(out, 'w')
	o.write('chr,window_start,window_end,rate\n')

	chr_lengths = { 'chr10': 20806668, 'chr11': 21403021, 'chr12': 21576510, 'chr13': 16962381,
                'chr14': 16419078, 'chr15': 14428146, 'chr16': 9909, 'chr17': 11648728,
                'chr18': 11201131, 'chr19': 11587733, 'chr1A': 73657157, 'chr1B': 1083483,
                'chr1': 118548696, 'chr20': 15652063, 'chr21': 5979137, 'chr22': 3370227,
                'chr23': 6196912, 'chr24': 8021379, 'chr25': 1275379, 'chr26': 4907541,
                'chr27': 4618897, 'chr28': 4963201, 'chr2': 156412533, 'chr3': 112617285,
                'chr4A': 20704505, 'chr4': 69780378, 'chr5': 62374962, 'chr6': 36305782,
                'chr7': 39844632, 'chr8': 27993427, 'chr9': 27241186, 'chrLG2': 109741,
                'chrLG5': 16416, 'chrLGE22': 883365, 'chrZ': 72861351 }

	chr_start = 0
	chr_end = chr_lengths[chr]

	for block_start in range(chr_start, chr_end, block_size):
                block_end = block_start + block_size
                if block_end > chr_end:
                        block_end = chr_end

                tmp = d[d.right_snp >= block_start]
                tmp = tmp[tmp.left_snp <= block_end]
                
                chunk_rho = {'block': {'bp': 0, 'rho': 0}}
                chunks = {'block': [int(block_start), int(block_end)]}
		
                for start, end, rate in izip(tmp.left_snp, tmp.right_snp, tmp.meanrho):
                        start = int(start)
                        end = int(end)

                        for chunk in chunks:
                                rate_mult = 0
                                diff = 0
                                # contained within
                                if start >= chunks[chunk][0] and end <= chunks[chunk][1]:
                                        rate_mult = (end - start) * rate
                                        diff = end - start
                                # hanging off to the left
                                elif start < chunks[chunk][0] and (end <= chunks[chunk][1] and end >= chunks[chunk][0]):
                                        rate_mult = (end - chunks[chunk][0]) * rate
                                        diff = end - chunks[chunk][0]
                           	# hanging off to the right
			   	elif (start >= chunks[chunk][0] and start <= chunks[chunk][1]) and end > chunks[chunk][1]:
                                        rate_mult = (chunks[chunk][1] - start) * rate
                                        diff = (chunks[chunk][1] - start)
				# hanging off on both sides
				elif start <= chunks[chunk][0] and end >= chunks[chunk][1]:
					diff = (chunks[chunk][1] - chunks[chunk][0])
					rate_mult = diff * rate
                                chunk_rho[chunk]['rho'] += rate_mult
                                chunk_rho[chunk]['bp'] += diff

		block_rate = 'NA'

		if chunk_rho['block']['bp'] > 0:
                	block_rate = '%.5f' % (chunk_rho['block']['rho'] / chunk_rho['block']['bp'])

                o.write('%s,%s,%s,%s\n' % (chr, block_start, block_end, block_rate))
	o.close()
	return
	
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--chr", help="chromosome for which to run analysis")
	parser.add_argument("--window", help="size in bp for window to smooth values")
	parser.add_argument("--sp", help="species [ZF|LTF] for which to run")

	args = parser.parse_args()
	chr = args.chr
	block_size = int(args.window)
	sp = args.sp
	
	dir = '/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/' % (sp)
	out = '%smaps/%s.window%s.bpen100.rm.txt' % (dir, chr, block_size)
	file = '%smaps/%s_recombination_bpen100.rm.txt' % (dir, chr)

	find_hotspots(file, out, chr, block_size) 

if __name__ == "__main__":
    main()

