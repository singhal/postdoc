import glob
import pandas as pd
import re

files = glob.glob('/mnt/gluster/home/sonal.singhal1/simulations/FDR/maps/*hotspots*')
counts = {}

for file in files:
	d = pd.read_csv(file)
	d.rename(columns={'diff':'diffrate'}, inplace=True)
	rho = re.search('rho([0-9|\.]+)', file).group(1)
	switch = re.search('switch([0-9|\.]+)', file).group(1)
	if rho not in counts:
		counts[rho] = {}
	if switch not in counts[rho]:
		counts[rho][switch] = {'Mb': 0, 'false_pos': 0}
	counts[rho][switch]['false_pos'] += d[d.diffrate > 10].shape[0]
	counts[rho][switch]['Mb'] += 1

print 'rho,switch,false_pos/Mb'
for rho in counts:
	for switch in counts[rho]:
		print '%s,%s,%s' % (rho, switch, counts[rho][switch]['false_pos'] / float(counts[rho][switch]['Mb']))
