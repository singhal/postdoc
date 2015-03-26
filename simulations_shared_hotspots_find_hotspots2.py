import re
import glob
import pandas as pd
import numpy as np
from scipy.stats.stats import pearsonr

block = 2000
max_dist = 3000
heat = 10
zf_maps = glob.glob('/mnt/gluster/home/sonal.singhal1/simulations/shared/ZF/maps/*hotspot*')

def get_hotspots(map):
	d = pd.read_csv(map)
	d = d.rename(columns={'diff': 'rho_lambda'})

	d = d[d.rho_lambda >= heat]

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

	hotspots = {}
	for ix, group in groups:
		start = group.block_start.min()
		length = group.block_end.max() - start
		hotspots[start] = {	'flank_rate': group.flank_rate.mean(), 'block_rate': group.block_rate.mean(),
					'length': length, 
					'lambda_heat': (group.block_rate.mean() / group.flank_rate.mean())}
	
	return hotspots

print 'background_rho,hotspot_lambda,sim_num,num_ZF_spots,num_LTF_spots,avg_dist,num_matched,corr_heats,corr_back_rhos'
for zf_map in zf_maps:
	rho = float(re.search('rho([0-9|\.]+)', zf_map).group(1))
	hot_heat = int(re.search('diff(\d+)', zf_map).group(1))
	sim_num = int(re.search('_(\d)_5_hot', zf_map).group(1))

	ltf_map = '/mnt/gluster/home/sonal.singhal1/simulations/shared/LTF/maps/recombination_rho%s_diff%s_%s_5_hotspots_blocksize2000_flanksize40000.txt' % (rho/ 2.0, hot_heat, sim_num)

	zf_hot = get_hotspots(zf_map)
	ltf_hot = get_hotspots(ltf_map)
	
	match = 0
	heats = {'ZF': [], 'LTF': []}
	rhos = {'ZF': [], 'LTF': []}
	starts = {'ZF': [], 'LTF': []}
	for start in zf_hot:
		if len(ltf_hot) > 0:
			nearest_matches = [abs(x - start) for x in sorted(ltf_hot.keys())]
			if np.min(nearest_matches) <= max_dist:
				match += 1
				ltf_match = sorted(ltf_hot.keys())[nearest_matches.index(np.min(nearest_matches))]
				heats['ZF'].append(zf_hot[start]['lambda_heat'])
				heats['LTF'].append(ltf_hot[ltf_match]['lambda_heat'])
				rhos['ZF'].append(zf_hot[start]['flank_rate'])
				rhos['LTF'].append(ltf_hot[ltf_match]['flank_rate'])
				starts['ZF'].append(start)
				starts['LTF'].append(ltf_match)

	if len(ltf_hot) > 0:
		avg_dist = np.mean([abs(x - y) for x, y in zip(starts['ZF'], starts['LTF'])])
		print '%s,%s,%s,%s,%s,%.0f,%s,%.3f,%.3f' % (	rho, hot_heat, sim_num, len(zf_hot), len(ltf_hot), avg_dist, match, 
							pearsonr(heats['ZF'], heats['LTF'])[0],
							pearsonr(rhos['ZF'], rhos['LTF'])[0] )
	else:
		print '%s,%s,%s,%s,0,NA,0,NA,NA' % (rho, hot_heat, sim_num, len(zf_hot))
