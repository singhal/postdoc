import pandas as pd
import numpy as np
import glob
import re

heat = 5
block = 2000

def get_hotspots(map):
        d = pd.read_csv(map)
        d = d.rename(columns={'diff': 'rho_lambda'})

        d = d[d.rho_lambda >= heat]

        d['dist'] = d.block_start.diff()
        groups = []
        group = -1
        for diff in d.dist:
                if diff > (block / 2.0) or pd.isnull(diff):
                        group += 1
                        groups.append(group)
                else:
                        groups.append(group)
        d['group'] = groups
        groups = d.groupby('group')

        hotspots = {}
        for ix, group in groups:
                start = group.block_start.min()
                length = group.block_end.max() - start
		mid = start + int(length / 2.0)
                hotspots[mid] = {     'flank_rate': group.flank_rate.mean(), 'block_rate': group.block_rate.mean(),
                                        'length': length, 
                                        'lambda_heat': (group.block_rate.mean() / group.flank_rate.mean())}
        
        return hotspots

o = open('/mnt/gluster/home/sonal.singhal1/simulations/shared/seqldhot/ZF_LTF.putative_hotspots.lambda5.csv', 'w')
files = glob.glob('/mnt/gluster/home/sonal.singhal1/simulations/shared/*/maps/*hotspots_blocksize2000_flanksize40000.txt')
o.write('file,midpoint,flank_rate,block_rate,length,lambda_heat\n')
keys = ['flank_rate', 'block_rate', 'length', 'lambda_heat']
rhos = {'0.0005': '0.001', '0.001': '0.002', '0.005': '0.01', '0.05': '0.1', '0.25': '0.5', '0.4': '0.8'}
for file in files:
	chr = re.search('(recombination_rho.*_5)_hotspots', file).group(1)
	if re.search('LTF', file):
		old_rho = re.search('rho([\d|\-|\.]+)', chr).group(1)
		new_rho = rhos[old_rho]
		chr = chr.replace(old_rho, new_rho)
	hotspots = get_hotspots(file)
	for mid in hotspots:
		o.write('%s,%s,%s\n' % (chr, mid, ','.join([str(hotspots[mid][key]) for key in keys])))
o.close()
