import re
import numpy as np
import pandas as pd
import scipy.stats

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

d = pd.read_csv('/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/spot2kb_flank40kb.seqldhot_validate_hotspots.csv')

z = d[d.zlk >= 10]

z['heat_bin'] = [heat_vals[x] for x in np.digitize(z.zheat, heat_bins)]
z['rho_bin'] = [rho_vals[x] for x in np.digitize(z.zbackrho, rho_bins)]

z = z.groupby(['heat_bin', 'rho_bin'])
for (heat_bin, rho_bin), group in z:
	num_spots_tested = group.shape[0]
	group = group[group.llk >= 10]
	num_matched = group.shape[0]
	per_matched = num_matched / float(num_spots_tested)

	corr = scipy.stats.pearsonr(group.zheat, group.lheat)
	dist = np.mean([abs(x-y) for x, y in zip(group.zmid, group.lmid)])

	print '%s,%s,%s,%.2f,%s,%.2f' % (heat_bin, rho_bin, num_spots_tested, per_matched, dist, corr[0])
	
