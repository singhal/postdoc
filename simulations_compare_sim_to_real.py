import glob
import re
import pandas as pd

dir = '/mnt/gluster/home/sonal.singhal1/simulations/hotspot_simulations/sim_bpen/'
seq_size = 1e6

tack_on = 'bpen'
outfile = '%s%s.compare_sim_to_real.csv' % (dir, tack_on)
o = open(outfile, 'w')
o.write('theta,actualrho,bpen,estimatedrho\n')

map_files = glob.glob('%smaps/*.txt' % (dir))
map_files = [x for x in map_files if not re.search('hotspots', x)]

for map_file in map_files:
	print map_file
	hotspot_d = re.search('(rho[0-9|\.]+_\d+)_bpen', map_file).group(1)
	hotspot_file = '/mnt/gluster/home/sonal.singhal1/simulations/hotspot_simulations/sim_10_20_40_60_80_100/hotspots/hotspot_%s.txt' % hotspot_d
	
	# get in hotspots
	hotspot_f = open(hotspot_file, 'r')
	hotspots = {}
	for l in hotspot_f:
		d = re.split('\t', l)
		d = [float(x) for x in d]
		for i in range(int(d[0] * seq_size), int(d[1] * seq_size) + 1):
			hotspots[i] = 1
	hotspot_f.close()

	theta = 0.0073
	type = re.search('_(rho.*)\.txt', hotspot_file).group(1)
	if re.search('theta([\d|\.]+)', hotspot_file):
		theta = float(re.search('theta([\d|\.]+)', hotspot_file).group(1))
		type = re.search('_(theta[^\/]*)\.txt', hotspot_file).group(1)
	rho = float(re.search('rho([\d|\.|\-|e]+)', hotspot_file).group(1))
	
	bpen = int(re.search('bpen(\d+)', map_file).group(1))
	
	d = pd.read_csv(map_file, sep=" ", skiprows=3, header=None, names=['left_snp', 'right_snp', 'meanrho', 'p025', 'p975'])
	for ix, (pos, meanrho) in enumerate(zip(d.left_snp, d.meanrho)):
		if ix % 40 == 0:
			if pos not in hotspots:
				o.write('%s,%s,%s,%s\n' % (theta, rho, bpen, meanrho))
