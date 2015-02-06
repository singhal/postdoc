import pandas as pd
import glob
import re
from itertools import izip
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--sp", help="species for which to run analysis")
parser.add_argument("--heat", help="heat cutoff to use")
args = parser.parse_args()
sp = args.sp
heat = int(args.heat)

# distance between two putative hotspots
# if the distance between two putative hotspots is less than this, don't keep both
dist = 5000

chrs = ['chr1', 'chr1A', 'chr2', 'chr3', 'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
	'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chrZ']

out = '/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/%s.putative_hotspots.heat%s.out' % (sp, sp, heat)
o = open(out, 'w')
o.write('chr,species,spot_start,block,block_rate,flank,flank_rate,lambda_heat,length\n')
for chr in chrs:
	hotspots = {}
	for flank in [40000,20000]:
		for block in [2000,1000,500]:	
			infile = '/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/putative_hotspots/%s.putativehotspots.block%s_flank%s.out' % (sp, chr, block, flank)
			d = pd.read_csv(infile)
			d = d[ d.rate_ratio > heat ]

			d['dist'] = d.block_start.diff()
			groups = []
			group = -1
			for diff in d.dist:
				if diff > (block / 2.0) or pd.isnull(diff):
					group += 1
					groups.append(group)
				else:
					groups.append(group)
			d['group'] = groups
			groups = d.groupby('group')

			for ix, group in groups:
				start = group.block_start.min()
				length = group.block_end.max() - start
				keep = True
				if len(hotspots) > 0:
					if np.min([abs(start - x) for x in hotspots]) < dist:
						keep = False
				if keep:
					hotspots[start] = {'flank_rate': '%.4f' % group.flank_rate.mean(), 'block_rate': '%.4f' % group.block_rate.mean(), 
							   'block': block, 'flank': flank, 'length': length, 'lambda_heat': '%.4f' % (group.block_rate.mean() / group.flank_rate.mean())}
	if hotspots:
		starts = sorted(hotspots.keys())
		columns = sorted(hotspots[starts[0]].keys())	
		for start in starts:
			o.write('%s,%s,%s,' % (chr, sp, start))
			o.write('%s' % ','.join([str(hotspots[start][column]) for column in columns]))
			o.write('\n')
o.close()
	

