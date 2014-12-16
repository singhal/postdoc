import re
import pandas as pd
import os
from itertools import izip
import glob

heat_cutoff = 10

def parse_file(file, heat_cutoff):
	posterior = None
	mean_heat = None
	mean_025 = None
	mean_975 = None
	mean_rho = None
	
	if (os.path.isfile(file)):
		if (os.path.getsize(file) > 0):
			d = pd.read_csv(file, header=None, sep="\s+",names=['rho','start','end','heat'])					
		
			posterior = d[d.heat > 10].shape[0] / float(d.shape[0])
			mean_heat = d.heat.mean()
			mean_025 = d.heat.quantile(q=0.025)
			mean_975 = d.heat.quantile(q=0.975)
			mean_rho = d.rho.mean()

	return posterior, mean_heat, mean_025, mean_975, mean_rho

files = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/analysis/phase_hotspots/rephased/*_hotspot')

print 'chr,start,as_posterior,as_meanheat,as_mean025,as_mean975,as_meanrho,re_posterior,re_meanheat,re_mean025,re_mean975,re_meanrho'
for rephased_file in files:
	asphased_file = rephased_file.replace('rephased/', '')
	asphased_file = asphased_file.replace('rephased', 'asphased')

	start = re.search('_(\d+)', rephased_file).group(1)
	chr = re.search('(chr[A-Z|0-9]+)', rephased_file).group(1)

	aposterior, amean_heat, amean_025, amean_975, amean_rho = parse_file(asphased_file, heat_cutoff)
	rposterior, rmean_heat, rmean_025, rmean_975, rmean_rho = parse_file(rephased_file, heat_cutoff)
	
	print '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s' % (start, chr, aposterior, amean_heat, amean_025, amean_975, amean_rho, \
							rposterior, rmean_heat, rmean_025, rmean_975, rmean_rho)
