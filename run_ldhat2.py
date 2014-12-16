import re
import glob
import subprocess

dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/hotspot_checks/LDhat/'
res_files = glob.glob(dir + 'results/*rates*')

for res_file in res_files:
	prefix = re.search('(\S+)rates', res_file).group(1)
	loc = res_file.replace('rates.txt', '.locs')
	loc = loc.replace('results', 'sites_locs')

	subprocess.call("~/bin/LDhat_v2.2/stat -input %s -burn 100 -loc %s -prefix %s" % (res_file, loc, prefix), shell=True)
