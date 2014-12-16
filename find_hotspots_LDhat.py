import pandas as pd
import argparse
import numpy as np
from itertools import izip
import os
import glob

def find_hotspots(file, out, block_size, flank_size):
	d = pd.read_csv(file, sep='\s+')[1:]
	d.Loci = d.Loci * 1000
	d['left_snp'] = d.Loci
	d['right_snp'] = d.Loci[1:].tolist() + [max(d.Loci[1:].tolist())]
	o = open(out, 'w')
	o.write('block_start,block_end,flank_rate,block_rate,rate_ratio\n')

	chr_start = int(d.Loci.min())
	chr_end = int(d.Loci.max())

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
		
                for start, end, rate in zip(tmp.left_snp, tmp.right_snp, tmp.Median):
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

                o.write('%s,%s,%s,%s,%s\n' % (block_start, block_end, flank_rate, block_rate, diff))
	o.close()
	return
	
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--block", help="size in bp for block size to evaluate as putative hotspot")
	parser.add_argument("--flank", help="size in bp for flank size to evaluate on either side of hotspot")	

	args = parser.parse_args()
	block_size = int(args.block)
	flank_size = int(args.flank)
	
	dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/hotspot_checks/LDhat/results/'
	files = glob.glob(dir + '*res*')
	for file in files:
		print file
		out = file.replace('res.txt', 'hotspots.txt')
		if not os.path.isfile(out):
			find_hotspots(file, out, block_size, flank_size) 

if __name__ == "__main__":
    main()

