import glob
import re
import pandas as pd

seq_size = 1e6
dir = '/mnt/gluster/home/sonal.singhal1/simulations/hotspot_simulations/sim_10_20_40_60_80_100/'
hotspot_files = glob.glob('%shotspots/*' % dir)

tack_on = '_blocksize2000_flanksize40000'
outfile = '%shotspot_detection%s.csv' % (dir, tack_on)
out = open(outfile, 'w')
out.write('real_rho,theta,real_length,bpen,actual_rate_diff,estimated_rate_diff,false_pos\n')

outfile1 = '%shotspot_power%s.csv' % (dir, tack_on)
out1 = open(outfile1, 'w')
out1.write('bpen,theta,rho,cutoff,power\n')

abs_cut_off = 10
cut_off = {10: 100, 20: 100, 40: 100, 60: 100, 80: 100, 100: 100}
# cut_off = {10: 48, 50: 48, 100: 48}
found_hotspots = {}

for hotspot_file in hotspot_files:
	# get in hotspots
	hotspot_f = open(hotspot_file, 'r')
	hotspots = {}
	for l in hotspot_f:
		d = re.split('\t', l)
		d = [float(x) for x in d]
		if d[2] not in hotspots:
			hotspots[d[2]] = []
		hotspots[d[2]].append([int(d[0] * seq_size), int(d[1] * seq_size)])
	hotspot_f.close()

	theta = 0.0073
	type = re.search('_(rho.*)\.txt', hotspot_file).group(1)
	if re.search('theta([\d|\.]+)', hotspot_file):
		theta = float(re.search('theta([\d|\.]+)', hotspot_file).group(1))
		type = re.search('_(theta.*)\.txt', hotspot_file).group(1)
	rho = float(re.search('rho([\d|\.|\-|e]+)', hotspot_file).group(1))
	maps = glob.glob('%smaps/*%s_*hotspots*%s*' % (dir, type, tack_on))

	if rho not in found_hotspots:
		found_hotspots[rho] = {}
	if theta not in found_hotspots[rho]:
		found_hotspots[rho][theta] = []

	for map in maps:
		bpen = int(re.search('(\d+)_hotspot', map).group(1))
		d = pd.read_csv(map)
		d = d.rename(columns = {'diff':'diff_rate'})
		# want to evaluate false positive, false negative, rate

		for spot in hotspots:
			num_false = 0
			# false positive
			poss_hotspots = d[d.diff_rate > spot]
			for hot_start, hot_rate in zip(poss_hotspots.block_start, poss_hotspots.diff_rate):
				false_positive = True
				for i in hotspots:
					for start, end in hotspots[i]:
						if abs(start - hot_start) < 1000:
							false_positive = False
				if false_positive:
					num_false += 1
			
			# false negative & rate
			for start, end in hotspots[spot]:
				diff_values = abs(d.block_start - start).tolist()
				ix = diff_values.index(min(diff_values))
				# this is probably the best match
				# note this is the older version of pandas indexing
				match = d[ix:ix+1]
				if match.diff_rate >= abs_cut_off:
					found_hotspots[rho][theta].append(spot)
				diff_rate = float(match.diff_rate) / spot
				out.write('%s,%s,%s,%s,%s,%.2f,%s\n' % (rho, theta, end - start, bpen, spot, diff_rate, num_false))
for rho in found_hotspots:
	for theta in found_hotspots[rho]:
		for cut in cut_off:
			count = found_hotspots[rho][theta].count(cut)
			out1.write('%s,%s,%s,%s,%s\n' % (bpen, theta, rho, cut, count / float(cut_off[cut])))
out.close()
out1.close()
				
