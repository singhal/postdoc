import pandas as pd
import argparse
import numpy as np
from itertools import izip

def find_hotspots(file, out, chr, block_size, flank_size):
	d = pd.read_csv(file, sep=" ", skiprows=3, header=None, 
		names=['left_snp', 'right_snp', 'meanrho', 'p0.05', 'p0.95'])
	o = open(out, 'w')
	o.write('chr,block_start,block_end,flank_rate,block_rate,rate_ratio\n')

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

	for block_start in range(chr_start, chr_end, int(block_size / 2.0)):
                block_end = block_start + block_size
                if block_end > chr_end:
                        block_end = chr_end

                left_flank_start = chr_start
                if (block_start - flank_size) >= left_flank_start:
                        left_flank_start = block_start - flank_size
                tmp = d[d.right_snp >= left_flank_start]

                right_flank_end = chr_end
                if (block_end + flank_size) <= right_flank_end:
                        right_flank_end = block_end + flank_size
                tmp = tmp[tmp.left_snp <= right_flank_end]
                
                chunk_rho = {'left': {'bp': 0, 'rho': 0}, 'block': {'bp': 0, 'rho': 0}, 'right': {'bp': 0, 'rho': 0}}
                chunks = {'left': [int(left_flank_start), int(block_start)], 'block': [int(block_start), int(block_end)], 'right': [int(block_end), int(right_flank_end)]}
		
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


		flank_rate = 'NA'
		block_rate = 'NA'
		diff = 'NA'
		if (chunk_rho['left']['bp'] + chunk_rho['right']['bp']) > 0:
			flank_rate = '%.5f' % ((chunk_rho['left']['rho'] + chunk_rho['right']['rho']) / float(chunk_rho['left']['bp'] + chunk_rho['right']['bp']))
		if chunk_rho['block']['bp'] > 0:
                	block_rate = '%.5f' % (chunk_rho['block']['rho'] / chunk_rho['block']['bp'])
		if block_rate != 'NA' and flank_rate != 'NA':	
			if float(flank_rate) > 0:
				diff = '%.5f' % (float(block_rate) / float(flank_rate))       

                o.write('%s,%s,%s,%s,%s,%s\n' % (chr, block_start, block_end, flank_rate, block_rate, diff))
	o.close()
	return
	
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--chr", help="chromosome for which to run analysis")
	parser.add_argument("--block", help="size in bp for block size to evaluate as putative hotspot")
	parser.add_argument("--flank", help="size in bp for flank size to evaluate on either side of hotspot")	
	parser.add_argument("--sp", help="species [ZF|LTF] for which to run")

	args = parser.parse_args()
	chr = args.chr
	block_size = int(args.block)
	flank_size = int(args.flank)
	sp = args.sp
	
	dir = '/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/with_fam/' % (sp)
	out = '%sputative_hotspots/%s.putativehotspots.block%s_flank%s.out' % (dir, chr, block_size, flank_size)
	file = '%smaps/%s_recombination_bpen5.txt' % (dir, chr)

	find_hotspots(file, out, chr, block_size, flank_size) 

if __name__ == "__main__":
    main()

