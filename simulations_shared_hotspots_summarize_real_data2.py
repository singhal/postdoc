import pandas as pd
import numpy as np
import scipy.stats

z = pd.read_csv('/mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/ZF.putative_hotspots.heat5.out')
z = z[z.flank == 40000]
z = z[z.block == 2000]

def get_vals(bins):
        vals = {}
        for ix, i in enumerate(bins[:-1]):
                vals[ix + 1] = '%s-%s' % (bins[ix], bins[ix + 1])
        vals[0] = '<%s' % bins[0]
        vals[np.max(vals.keys())+1] = '>%s' % (bins[-1])
        return vals

tranches = {}
heat_bins = [10, 15, 20] 
rho_bins = [0.002, 0.01, 0.1, 0.5]

heat_vals = get_vals(heat_bins)
rho_vals = get_vals(rho_bins)

z['heat_bin'] = [heat_vals[x] for x in np.digitize(z.lambda_heat, heat_bins)]
z['rho_bin'] = [rho_vals[x] for x in np.digitize(z.flank_rate, rho_bins)]

l = pd.read_csv('/mnt/gluster/home/sonal.singhal1/LTF/analysis/LDhelmet/LTF.putative_hotspots.heat5.out')

tranches = {}
for i in heat_vals.values():
        tranches[i] = {}
        for j in rho_vals.values():
                tranches[i][j] = {'num_spots': 0, 'num_matched': 0, 'dists': [], 'heat_ZF': [], 'heat_LTF': []}


z = z.groupby(['chr'])
for chr, groupz in z:
	groupl = l[l.chr == chr]
	ltf_starts = groupl.spot_start.tolist()

	for heat_bin, rho_bin, zstart, zheat in zip(groupz.heat_bin, groupz.rho_bin, groupz.spot_start, groupz.lambda_heat):
		tranches[heat_bin][rho_bin]['num_spots'] += 1
		
		dists = [abs(x - zstart) for x in ltf_starts]
		min_dist = np.min(dists)
		if min_dist < 3000:
			tranches[heat_bin][rho_bin]['num_matched'] += 1
			match = ltf_starts[dists.index(min_dist)]

			tranches[heat_bin][rho_bin]['dists'].append(min_dist)
			tranches[heat_bin][rho_bin]['heat_ZF'].append(zheat)
			tranches[heat_bin][rho_bin]['heat_LTF'].append(float(groupl[groupl.spot_start == match].lambda_heat))

for heat_bin in tranches:
        for rho_bin in tranches[heat_bin]:
                if tranches[heat_bin][rho_bin]['num_spots'] > 0:
                        per_matches = tranches[heat_bin][rho_bin]['num_matched'] / float(tranches[heat_bin][rho_bin]['num_spots'])
                        avg_dist = np.mean(tranches[heat_bin][rho_bin]['dists'])
                        corr_heat = scipy.stats.pearsonr(tranches[heat_bin][rho_bin]['heat_ZF'], tranches[heat_bin][rho_bin]['heat_LTF'])
                        print '%s,%s,%s,%.2f,%.0f,%.2f' % (heat_bin, rho_bin, tranches[heat_bin][rho_bin]['num_spots'], per_matches, avg_dist, corr_heat[0])


