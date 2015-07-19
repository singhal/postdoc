import glob
import pandas as pd
import re

files = glob.glob('/mnt/gluster/home/sonal.singhal1/simulations/FDR/maps/*hotspots*')
counts = {}
rho_lambda = 5

def removeCloseItems(items, itemDistance):
    if items:
        lastOutput = items[0]
        yield items[0]
        for currentItem in items[1:]:
            if ((currentItem - lastOutput) > itemDistance):
                lastOutput = currentItem
                yield currentItem

for file in files:
	d = pd.read_csv(file)
	d.rename(columns={'diff':'diffrate'}, inplace=True)
	rho = re.search('rho([0-9|\.]+)', file).group(1)
	switch = re.search('switch([0-9|\.]+)', file).group(1)
	if rho not in counts:
		counts[rho] = {}
	if switch not in counts[rho]:
		counts[rho][switch] = {'Mb': 0, 'false_pos': 0}
	starts = d[d.diffrate >= rho_lambda].block_start.tolist()
	if len(starts) > 1:
		starts = list(removeCloseItems(starts, 2000))
	
	counts[rho][switch]['false_pos'] += len(starts)
	counts[rho][switch]['Mb'] += 1


print 'rho,switch,false_pos/Mb'
for rho in counts:
	for switch in counts[rho]:
		print '%s,%s,%s' % (rho, switch, counts[rho][switch]['false_pos'] / float(counts[rho][switch]['Mb']))
