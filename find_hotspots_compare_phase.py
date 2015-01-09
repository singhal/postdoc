import pandas as pd
import glob
import re
from itertools import izip

# want hotspots with at least this much heat
heat = 10
# if two species have hotspots within this overlap, we will say they share the hotspot
overlap = 5000

hotspots = {}
sp = 'ZF'

for block in [2000]:
	hotspots[block] = {}
	for flank in [40000]:
		hotspots[block][flank] = {}
		for subdir in ['putative_hotspots', 'old_putative_hotspots']:
			dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/%s/' % subdir
			files = glob.glob('%s*block%s_flank%s*' % (dir, block, flank))
			for file in files:
				chr = re.search('(chr[A-Z|a-z|0-9]+)', file).group(1)
				if chr not in hotspots[block][flank]:
					hotspots[block][flank][chr] = {}
				if subdir not in hotspots[block][flank][chr]:
					hotspots[block][flank][chr][subdir] = {}
			
				d = pd.read_csv(file)
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
					hotspots[block][flank][chr][subdir][start] = {'flank': group.flank_rate.mean(), 'block': group.block_rate.mean(), 'length': length}


def print_hotspots(hotspots, sd1, sd2):
	for block in hotspots:
		for flank in hotspots[block]:
	        	for chr in hotspots[block][flank]:
	                	if sd1 in hotspots[block][flank][chr]:
	                        	for spot in hotspots[block][flank][chr][sd1].keys():                                
						if sd2 in hotspots[block][flank][chr]:
	                                        	dist = filter(lambda x: abs(spot - x) < overlap, hotspots[block][flank][chr][sd2].keys())
	                                        	if dist:                                                
								o.write('%s,%s,%s,%s,%s,%s,%.4f,%.4f,%s,%s,%s,%.4f,%.4f\n' % (sd1, block, flank, chr, spot, 
									hotspots[block][flank][chr][sd1][spot]['length'],
									hotspots[block][flank][chr][sd1][spot]['block'],
									hotspots[block][flank][chr][sd1][spot]['flank'],
									'YES', dist[0], hotspots[block][flank][chr][sd2][dist[0]]['length'],
									hotspots[block][flank][chr][sd2][dist[0]]['block'],
									hotspots[block][flank][chr][sd2][dist[0]]['flank']))
							else:   
        	                                		o.write('%s,%s,%s,%s,%s,%s,%.4f,%.4f,%s,%s,%s,%s,%s\n' % (sd1, block, flank, chr, spot, 
										hotspots[block][flank][chr][sd1][spot]['length'],
                                                                		hotspots[block][flank][chr][sd1][spot]['block'],
                                                                		hotspots[block][flank][chr][sd1][spot]['flank'],  'NO', 'NA', 'NA',
										'NA','NA'))
        	                        	else:
							o.write('%s,%s,%s,%s,%s,%s,%.4f,%.4f,%s,%s,%s,%s,%s\n' % (sd1, block, flank, chr, spot,
									hotspots[block][flank][chr][sd1][spot]['length'],
									hotspots[block][flank][chr][sd1][spot]['block'],
									hotspots[block][flank][chr][sd1][spot]['flank'],  'NO', 'NA', 'NA', 'NA', 'NA'))
	return

o = open('ZF.putative_hotspots.phase_compare.csv', 'w')
o.write('dir,spot_size,flank_size,chr,spot_start,length,block_rate,flank_rate,match_in_other_species,match_spot_start,match_length,match_block_rate,match_flank_rate\n')
print_hotspots(hotspots, 'putative_hotspots', 'old_putative_hotspots')
print_hotspots(hotspots, 'old_putative_hotspots', 'putative_hotspots')
o.close()

