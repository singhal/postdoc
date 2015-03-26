import pandas as pd
import numpy as np
from itertools import izip

zf_dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/'
ltf_dir = '/mnt/gluster/home/sonal.singhal1/LTF/analysis/hotspots/'

file = '%sspot2kb_flank40kb.seqldhot_validate_hotspots.csv' % (zf_dir)

d = pd.read_csv(file)
d = d[d.llk >= 10]
d['shared'] = ['no'] * d.shape[0]
d['percentile'] = ['NA'] * d.shape[0]

percents = [0,25,50,75,100]
d.ix[ (d.zlk >= 10), 'shared' ] = 'yes'
percentiles = np.percentile(d.lheat.tolist(), percents)

for (low, high, lowvalue, highvalue) in izip(percents, percents[1:], percentiles, percentiles[1:]):
	d.ix[ (d.lheat >= lowvalue) & (d.lheat < highvalue), 'percentile'] = '%s_%s' % (low, high)
d.ix[(d.percentile == 'NA'), 'percentile'] = '%s_%s' % (percents[-2], percents[-1])

d = d.drop(['zheat', 'zlk', 'zlength', 'zstart'], axis=1)
d = d.rename(columns={'lheat': 'lambda', 'lstart': 'start', 'llength': 'length', 'llk': 'likelihood'})

d.to_csv('%sLTF.seqldhot_hotspots.heat5.csv' % ltf_dir, index=False)
	

