import pandas as pd
import re
import numpy as np

dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/'
lr = '%shotspots.LR_seqldhot.ZF_LTF.csv' % dir
sp = '%sspot2kb_flank40kb.seqldhot_validate_hotspots.csv' % dir

lr = pd.read_csv(lr)
sp = pd.read_csv(sp)

lr['zmid'] = np.nan
lr['lmid'] = np.nan

for ix, (chr, start) in enumerate(zip(lr.chr, lr.start)):
	tmp = sp[sp.chr == chr]
	dist1 = [abs(start - mid) for mid in tmp.zmid]
	dist2 = [abs(start - mid) for mid in tmp.lmid]

	min1 = np.nanmin(dist1)
	min2 = np.nanmin(dist2)

	if min1 < 5000 and min2 < 5000:
		if min1 < min2:
			pick = tmp.index[dist1.index(min1)]
		else:
			pick = tmp.index[dist2.index(min2)]
		lr.ix[ix, 'zmid'] = tmp.ix[pick, 'zmid']
		lr.ix[ix, 'lmid'] = tmp.ix[pick, 'lmid']
	elif min1 < 5000:
		pick = tmp.index[dist1.index(min1)]
		lr.ix[ix, 'zmid'] = tmp.ix[pick, 'zmid']
                lr.ix[ix, 'lmid'] = tmp.ix[pick, 'lmid']
	elif min2 < 5000:
		pick = tmp.index[dist2.index(min2)]
                lr.ix[ix, 'zmid'] = tmp.ix[pick, 'zmid']
                lr.ix[ix, 'lmid'] = tmp.ix[pick, 'lmid']
	else:
		#print lr.ix[ix]
		pass

lr.to_csv('%shotspots.LR_seqldhot.ZF_LTF_mids.csv' % dir, index=False)
