import pandas as pd
import argparse
import numpy as np
from itertools import izip

def find_hotspots(dir, out, chr, block_size, flank_size, fold_diff):
	file = '%s%s.bpen10.txt' % (dir, chr)
	d = pd.read_csv(file, sep=" ", skiprows=3, header=None, 
		names=['left_snp', 'right_snp', 'meanrho', 'p0.025', 'p0.5', 'p0.975'])
	o = open(out, 'w')

	chr_start = min(d.left_snp)
	chr_end = max(d.right_snp)

	for block_start in range(chr_start, chr_end, block_size):
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
    
	    rec_rate = {}
	    for pos_start, pos_end, rate in izip(tmp.left_snp, tmp.right_snp, tmp.meanrho):
	        for pos in range(pos_start, pos_end):
	            rec_rate[pos] = rate
               
	    left_flank = [rec_rate[i] for i in range(left_flank_start, block_start)]
	    right_flank = [rec_rate[i] for i in range(block_end, right_flank_end)]
	    flank = left_flank + right_flank
    
	    block = [rec_rate[i] for i in range(block_start, block_end)]
    
	    if np.mean(block) >= np.mean(flank) * fold_diff:
	        o.write('%s,%s,%s,%s\n' % (block_start, block_end, np.mean(flank), np.mean(block)))
	o.close()
	return
	
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--chr", help="chromosome for which to run analysis")
	parser.add_argument("--block", help="size in bp for block size to evaluate as putative hotspot")
	parser.add_argument("--flank", help="size in bp for flank size to evaluate on either side of hotspot")	
	parser.add_argument("--fold_diff", help="fold difference for block rho to be higher than flank rho to call it a hotspot")

	args = parser.parse_args()
	chr = args.chr
	block_size = int(args.block)
	flank_size = int(args.flank)
	fold_diff = int(args.fold_diff)

	dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/drafts_Aug30/'
	out = '%s%s.putativehotspots.block%s_flank%s_diff%s.out' % (dir, chr, block_size, flank_size, fold_diff)

	find_hotspots(dir, out, chr, block_size, flank_size, fold_diff) 

if __name__ == "__main__":
    main()

