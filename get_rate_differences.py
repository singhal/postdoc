import pandas as pd

zf_dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/maps/'
ltf_dir = '/mnt/gluster/home/sonal.singhal1/LTF/analysis/LDhelmet/maps/'
cm_chrs = {'chr1': 55.126818712893566,
 'chr10': 74.116041662420571,
 'chr11': 43.960887850502928,
 'chr12': 27.339799229340038,
 'chr13': 17.821682598105092,
 'chr14': 41.038370882409204,
 'chr15': 50.805605798260906,
 'chr1A': 89.910986819301911,
 'chr2': 56.244681787471364,
 'chr3': 72.510492654788251,
 'chr4': 32.107282090789425,
 'chr4A': 35.274541314342983,
 'chr5': 88.402160152275684,
 'chr6': 71.877717562562324,
 'chr7': 59.515430582507577,
 'chr8': 41.301712170604233,
 'chr9': 58.280188040977514}

chrL = []
startL = []
stopL = []
zf_rhoL = []
ltf_rhoL = []
zf_cmL = []
ltf_cmL = []
rate_ratioL = []
rate_ratio_normalizedL = []

median_diff = 0.44

for chr in cm_chrs:
	zf_rho = pd.read_csv('%s%s.window100000.bpen100.txt' % (zf_dir, chr))
	# convert rho to cM / Mb
	ratio = (zf_rho.rho * (zf_rho.window_end - zf_rho.window_start)).cumsum().max() / cm_chrs[ chr ]
	zf_rho[ 'cmrate' ] = zf_rho.rho * (zf_rho.window_end - zf_rho.window_start) / ratio
	
	ltf_rho = pd.read_csv('%s%s.window100000.bpen100.txt' % (ltf_dir, chr))
	# convert rho to cM / Mb
	ratio = (ltf_rho.rho * (ltf_rho.window_end - ltf_rho.window_start)).cumsum().max() / cm_chrs[ chr ]
	ltf_rho[ 'cmrate' ] = ltf_rho.rho * (ltf_rho.window_end - ltf_rho.window_start) / ratio 

	# find relative difference
	zf_rho['rho_diff_unstandard'] = zf_rho.cmrate / ltf_rho.cmrate
	diff_standard = []
	# standardize with respect to mean
	for diff in zf_rho.rho_diff_unstandard:
		if diff < median_diff:
			diff_standard.append((-1 * median_diff) / diff)
		else:
			diff_standard.append( diff / median_diff)

	
	chrL += [chr] * len(diff_standard)
	startL += zf_rho.window_start.tolist()
	stopL += zf_rho.window_end.tolist()
	zf_rhoL += zf_rho.rho.tolist()
	ltf_rhoL += ltf_rho.rho.tolist()
	zf_cmL += zf_rho.cmrate.tolist()
	ltf_cmL += ltf_rho.cmrate.tolist()
	rate_ratioL += zf_rho.rho_diff_unstandard.tolist()
	rate_ratio_normalizedL += diff_standard

d = pd.DataFrame({'chr': chrL, 'window_start': startL, 'window_end': stopL, 'ZF_rho': zf_rhoL, 'LTF_rho': ltf_rhoL, 'ZF_cM_Mb': zf_cmL, 'LTF_cM_Mb': ltf_cmL, 'cM_MB_rate_ratio': rate_ratioL, 'rate_ratio_normalized': rate_ratio_normalizedL})

d.to_csv('/mnt/gluster/home/sonal.singhal1/for_ellen/ZF_LTF.recombination_rate_diffs.window100000.csv', index=False)
