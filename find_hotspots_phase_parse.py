import re
import pandas as pd
import os
from itertools import izip

putative_hotspots = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/LTF_ZF.putative_hotspots.csv'
chrs = ['chr1', 'chr1A', 'chr2', 'chr3', 'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', \
                'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chrZ']
out = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/phase_hotspots/spot2kb_flank40kb.phase_validate_hotspots.csv'

heat_cutoff = 10

# get the most conservative set of hotspots 
d = pd.read_csv(putative_hotspots)
d = d[(d.spot_size == 2000) & (d.flank_size == 40000)]
d = d[d.chr.isin(chrs)]
d['rho_lambda'] = d.block_rate / d.flank_rate

hotspots = d.to_dict()
indices = hotspots['chr'].keys()

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

o = open(out, 'w')
o.write(','.join(sorted(hotspots.keys())) + ',zposterior,zmean_heat,zmean_025,zmean_975,zmean_rho,lposterior,lmean_heat,lmean_025,lmean_975,lmean_rho\n')	
for index in indices:
	start = hotspots['spot_start'][index]
	chr = hotspots['chr'][index]
	rho_lambda = hotspots['rho_lambda'][index]
	block_rate = hotspots['block_rate'][index]

	out_zf = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/phase_hotspots/putative_hotspot_%s_%s.PHASE.asphased.out_hotspot' % (chr, start)
	out_ltf = '/mnt/gluster/home/sonal.singhal1/LTF/analysis/phase_hotspots/putative_hotspot_%s_%s.PHASE.asphased.out_hotspot' % (chr, start)	

	zposterior, zmean_heat, zmean_025, zmean_975, zmean_rho = parse_file(out_zf, heat_cutoff)
	lposterior, lmean_heat, lmean_025, lmean_975, lmean_rho = parse_file(out_ltf, heat_cutoff)
	o.write(','.join([str(hotspots[key][index]) for key in sorted(hotspots)]))
	o.write(',%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (zposterior, zmean_heat, zmean_025, zmean_975, zmean_rho, \
							lposterior, lmean_heat, lmean_025, lmean_975, lmean_rho))
o.close()
