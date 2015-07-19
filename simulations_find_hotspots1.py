import pandas as pd
import argparse
import re
import numpy as np
from itertools import izip
import glob
import os

def find_hotspots(file, block_size, flank_size):
	d = pd.read_csv(file, sep=" ", skiprows=3, header=None, names=['left_snp', 'right_snp', 'meanrho', 'p025', 'p975'])
	outfile = file.replace('.txt', '_hotspots_blocksize%s_flanksize%s.txt' % (block_size, flank_size))
	out = open(outfile, 'w')
	out.write('block_start,block_end,flank_rate,block_rate,diff\n')

	chr_start = min(d.left_snp)
	chr_end = max(d.right_snp)

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

		flank_rate = (chunk_rho['left']['rho'] + chunk_rho['right']['rho']) / float(chunk_rho['left']['bp'] + chunk_rho['right']['bp'])
		block_rate = chunk_rho['block']['rho'] / chunk_rho['block']['bp']
		diff = block_rate / flank_rate		

		out.write('%s,%s,%.3f,%.3f,%.1f\n' % (block_start, block_end, flank_rate, block_rate, diff))
    	out.close()
	return
	
def main():
	block_size = 2000
	flank_size = 40000

	files = glob.glob('/mnt/gluster/home/sonal.singhal1/simulations/shared/*/maps/*txt')
	files = filter(lambda x: not re.search('hotspots', x), files)
	for file in files:
		out = file.replace('.txt', '_hotspots_blocksize%s_flanksize%s.txt' % (block_size, flank_size))
		out1 = file.replace('.txt', '_hotspots.txt')
		if not os.path.isfile(out) and not os.path.isfile(out1):
			print file
			find_hotspots(file, block_size, flank_size) 

if __name__ == "__main__":
    main()

