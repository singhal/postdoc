import re
import pandas as pd
import glob
import os
import numpy as np
import scipy.stats

match_dist = 3000
block = 2000
heat = 5
zf_files = glob.glob('/mnt/gluster/home/sonal.singhal1/simulations/shared/ZF/maps/*hot*')

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
                mid = int((group.block_start.min() + group.block_end.max()) / 2.0)
                length = group.block_end.max() - group.block_start.min()
                hotspots[mid] = {     'flank_rate': group.flank_rate.mean(), 'block_rate': group.block_rate.mean(),
                                        'length': length, 
                                        'lambda_heat': (group.block_rate.mean() / group.flank_rate.mean())}
        
        return hotspots

def get_vals(bins):
	vals = {}
	for ix, i in enumerate(bins[:-1]):
        	vals[ix + 1] = '%s-%s' % (bins[ix], bins[ix + 1])
	vals[0] = '<%s' % bins[0]
	vals[np.max(vals.keys())+1] = '>%s' % (bins[-1])
	return vals

# organize data format
tranches = {}
heat_bins = [10, 15, 20] 
rho_bins = [0.005, 0.01, 0.1, 0.5]

heat_vals = get_vals(heat_bins)
rho_vals = get_vals(rho_bins)

for i in range(len(heat_bins)+2):
	tranches[i] = {}
	for j in range(len(rho_bins)+2):
		tranches[i][j] = {'num_spots': 0, 'num_matched': 0, 'dists': [], 'heat_ZF': [], 'heat_LTF': []} 

for zf_file in zf_files:
	ltf_file = zf_file.replace('ZF', 'LTF')
	zf_rho = float(re.search('rho([0-9|\.|-]+)', zf_file).group(1))
	ltf_rho = zf_rho / 2.0
	ltf_file = ltf_file.replace(str(zf_rho), str(ltf_rho))

	if os.path.isfile(zf_file):
		if os.path.isfile(ltf_file):
			zf_hot = get_hotspots(zf_file)
			ltf_hot = get_hotspots(ltf_file)
			ltf_starts = ltf_hot.keys()

			for spot in zf_hot:
				heat_bin = np.digitize([zf_hot[spot]['lambda_heat']], heat_bins)[0]
				rho_bin = np.digitize([zf_hot[spot]['flank_rate']], rho_bins)[0]

				tranches[heat_bin][rho_bin]['num_spots'] += 1
				
				dists = [abs(x - spot) for x in ltf_starts]
				if len(dists) > 0:
					if np.min(dists) < match_dist:
						tranches[heat_bin][rho_bin]['num_matched'] += 1	
						match = ltf_starts[dists.index(np.min(dists))]
						dist = abs(match - spot)
						tranches[heat_bin][rho_bin]['dists'].append(dist)
	
						tranches[heat_bin][rho_bin]['heat_ZF'].append(zf_hot[spot]['lambda_heat'])
						tranches[heat_bin][rho_bin]['heat_LTF'].append(ltf_hot[match]['lambda_heat'])

for heat_bin in tranches:
	for rho_bin in tranches[heat_bin]:
		if tranches[heat_bin][rho_bin]['num_spots'] > 0:
			per_matches = tranches[heat_bin][rho_bin]['num_matched'] / float(tranches[heat_bin][rho_bin]['num_spots'])
			avg_dist = np.mean(tranches[heat_bin][rho_bin]['dists'])
			corr_heat = scipy.stats.pearsonr(tranches[heat_bin][rho_bin]['heat_ZF'], tranches[heat_bin][rho_bin]['heat_LTF'])
			print '%s,%s,%s,%.2f,%.0f,%.2f' % (heat_vals[heat_bin], rho_vals[rho_bin], tranches[heat_bin][rho_bin]['num_spots'], per_matches, avg_dist, corr_heat[0])
