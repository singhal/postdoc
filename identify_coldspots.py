import re
import pandas as pd
import glob
import numpy as np
import argparse
from itertools import izip

parser = argparse.ArgumentParser()
parser.add_argument("--sp", help="species for which to run analysis")
args = parser.parse_args()
sp = args.sp

chrs = ['chr1', 'chr1A', 'chr2', 'chr3', 'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', \
                'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chrZ']

hot_dir = '/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/putative_hotspots/' % sp

putative_hotspots = glob.glob('/mnt/gluster/home/sonal.singhal1/*/analysis/LDhelmet/*heat5*out')
out = '/mnt/gluster/home/sonal.singhal1/%s/analysis/hotspots/%s.ldhelmet_unvalidated_coldspots.csv' % (sp, sp)
o = open(out, 'w')
o.write('chr,spot_start,spot_end,rho\n')

# get the most conservative set of hotspots
hotspots = {}
for putative_hotspot in putative_hotspots: 
        d = pd.read_csv(putative_hotspot)
        d = d[(d.block == 2000) & (d.flank == 40000)]
        for chr, start in zip(d.chr, d.spot_start):
                if chr not in hotspots:
                        hotspots[chr] = []
                if len(hotspots[chr]) > 5000:
                        min_dist = np.min([abs(x - start) for x in hotspots[chr]])
                        if min_dist > 0:
                                hotspots[chr].append(start)
                else:
                        hotspots[chr].append(start)

for chr in chrs:
	file = '%s%s.putativehotspots.block2000_flank40000.out' % (hot_dir, chr)
	d = pd.read_csv(file)

	d = d[d.rate_ratio > 0.9]
	d = d[d.rate_ratio < 1.1]
	d = d[d.flank_rate > 0.001]
	d = d[d.flank_rate < 0.1]
	d = d[np.isfinite(d.block_rate)]
	d = d[np.isfinite(d.flank_rate)]
	
	for start, end, rate in izip(d.block_start, d.block_end, d.flank_rate):
		min_dist = np.min([abs(x - start) for x in hotspots[chr]])
		if min_dist > 25000:
			o.write('%s,%s,%s,%s\n' % (chr, start, end, rate))
o.close() 

	
