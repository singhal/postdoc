import pandas as pd
import scipy.stats
import numpy as np

d = pd.read_csv('phasing_uncertainty.csv')
d = d[d.lk >= 10]
d['rel_start'] = d.putative_start - d.start

spots = d.groupby(['chr', 'putative_start'])

var_length = []
var_heat = []
var_start = []

for ix, spot in spots:
	var_length.append(scipy.stats.variation(spot.length))
	var_start.append(scipy.stats.variation(spot.start))
	var_heat.append(scipy.stats.variation(spot.heat))

print np.mean([x for x in var_start if np.isfinite(x)])
print np.mean([x for x in var_length if np.isfinite(x)])
print np.mean([x for x in var_heat if np.isfinite(x)])
