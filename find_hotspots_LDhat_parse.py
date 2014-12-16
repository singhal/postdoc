import glob
import re
import pandas as pd

files = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/hotspot_checks/LDhat/results/*hotspot*')

print 'chromosome,hotspot_start,ratio_inferred'
for file in files:
	d = pd.read_csv(file)
	hot_start = int(re.search('(\d+)hot', file).group(1)) 
	chr = re.search('(chr[\d|A-Z]+)', file).group(1)
	
	start  = hot_start - 5000
	end = hot_start + 5000
	d = d[d.block_start >= start]
	d = d[d.block_end <= end]

	print '%s,%s,%s' % (chr, hot_start, d.rate_ratio.max())
