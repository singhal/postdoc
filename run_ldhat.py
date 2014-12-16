import glob
import subprocess
import re
import os.path

dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/hotspot_checks/LDhat/'
files = glob.glob(dir + 'sites_locs/*locs')

for loc_file in files:
	site_file = loc_file.replace('.locs', '.sites')
	prefix = re.search('(chr.*_\d+)', site_file).group(1)

	expected_out = dir + 'results/' + prefix + 'rates.txt'

	if not os.path.isfile(expected_out):
		call = '/mnt/lustre/home/sonal.singhal1/bin/LDhat_v2.2/interval -seq %s -loc %s -lk /mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/hotspot_checks/LDhat/new_lk.txt -its 5000000 -samp 5000 -bpen 5.0 -prefix %s' % (site_file, loc_file, dir + 'results/' + prefix)
		print call
