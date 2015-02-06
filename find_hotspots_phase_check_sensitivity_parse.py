import re
import pandas as pd
import os
from itertools import izip
import glob

heat_cutoff = 10

def parse_file(file, heat_cutoff):
        posterior = 'NA'
        mean_heat = 'NA'
        mean_025 = 'NA'
        mean_975 = 'NA'
        mean_rho = 'NA'
        
        if (os.path.isfile(file)):
                if (os.path.getsize(file) > 0):
                        d = pd.read_csv(file, header=None, sep="\s+",names=['rho','start','end','heat'])                                        
                
                        posterior = d[d.heat > heat_cutoff].shape[0] / float(d.shape[0])
                        mean_heat = d.heat.mean()
                        mean_025 = d.heat.quantile(q=0.025)
                        mean_975 = d.heat.quantile(q=0.975)
                        mean_rho = d.rho.mean()

        return posterior, mean_heat, mean_025, mean_975, mean_rho

print 'chr,bp,match_or_not,posterior,mean_heat,mean_heat025,mean_heat975,mean_rho_inferred,rho_ratio_inferred_given'
files = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/hotspot_checks/phase_sensitivity/*hotspot')
for file in files:
	# get the rho that we had inferred from LDhelmet
	prior_file = file.replace('.asphased.out_hotspot', '.prior')
	f = open(prior_file, 'r')
	given_flank_rho = float(f.next().rstrip())

	chr = re.search('(chr[A-Z|0-9]+)', file).group(1)
	bp = re.search('(\d+)\.PHASE', file).group(1)

	posterior, mean_heat, mean_025, mean_975, mean_rho = parse_file(file, heat_cutoff)

	if posterior > 0.9 or mean_heat > heat_cutoff:
		match = True
	else:
		match = False
	rho_ratio = mean_rho / given_flank_rho
	print '%s,%s,%s,%s,%s,%s,%s,%s,%s' % (chr, bp, match, posterior, mean_heat, mean_025, mean_975, mean_rho, rho_ratio)

