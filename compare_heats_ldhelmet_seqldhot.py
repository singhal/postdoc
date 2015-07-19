import re
import pandas as pd
from itertools import izip
import numpy as np

d = pd.read_csv('/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/spot2kb_flank40kb.seqldhot_validate_hotspots.csv')
o = open('/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/correlation_btn_seqldhot_ldhelmet.csv', 'w')

z = d[d.zlk >= 10]
l = d[d.llk >= 10]

def get_rho_window(r, start, end):
	tmp = r[r.right_snp >= start]
	tmp = tmp[tmp.left_snp <= end]

	allrho = 0
	bp = 0	

	for block_start, block_end, rho in zip(tmp.left_snp, tmp.right_snp, tmp.meanrho):
		if block_start < start and block_end < end:
			bp += (block_end - start)
			allrho += rho * (block_end - start)
		elif block_start < start and block_end > end:
			bp += (end - start)
			allrho += rho * (end - start)
		elif block_start > start and block_end < end:
			bp += (block_end - block_start)
			allrho += rho * (block_end - block_start)
		elif block_start > start and block_end > end:
			bp += (end - block_start)
			allrho += rho * (end - block_start)

	if bp  > 0:
		return allrho / float(bp)
	else:
		return np.nan

def compare_heats(d, species, prefix):
	d = d.groupby('chr')
	for chr, g in d:
		r = '/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/maps/%s_recombination_bpen5.txt' % (species, chr)
		r = pd.read_csv(r, sep=" ", skiprows=3, header=None, 
					names=['left_snp', 'right_snp', 'meanrho', 'p0.05', 'p0.95'])
		
		for mid, length, seq_backrho, seq_heat in izip(g[prefix + 'mid'], g[prefix + 'length'], g[prefix + 'backrho'], g[prefix + 'heat']):
			hot_start = mid - int(length / 2.0)
			hot_end = mid + int(length / 2.0)

			start = hot_start - 40000
			end = hot_end + 40000

			hotspot_rho = get_rho_window(r, hot_start, hot_end)
			flank_left = get_rho_window(r, start, hot_start)
			flank_right = get_rho_window(r, hot_end, end)
		
			flanks = [flank_left, flank_right]
			ldh_backrho = np.mean([x for x in flanks if np.isfinite(x)])

			o.write('%s,%s,%s,%s,%s,%s,%s\n' % (chr, mid, species, seq_backrho, seq_heat, ldh_backrho, hotspot_rho))

o.write('chr,spot_midpoint,species,seqldhot_backrho,seqldhot_heat,ldh_backrho,ldh_hotspotrho\n')
compare_heats(z, 'ZF', 'z')
compare_heats(l, 'LTF', 'l')
o.close()
