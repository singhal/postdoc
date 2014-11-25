import pandas as pd
import glob
import re

files = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/maps/*.window10000.bpen100.rm.txt')

print 'chr,median_rate,mean_rate'
for file in files:
	chr = re.search('(chr[A-Z|0-9]+)', file).group(1)
	d = pd.read_csv(file)
	print '%s,%s,%s' % (chr, d.rate.median(), d.rate.mean())
