import glob
import subprocess
import os
import re

files = glob.glob('/mnt/gluster/home/sonal.singhal1/simulations/shared/*/maps/*')
files = filter(lambda x: not re.search('txt', x), files)
for file in files:
	out = '%s.txt' % file
	if not os.path.isfile(out):
		# print file
		subprocess.call('~/bin/LDhelmet_v1.6/ldhelmet post_to_text -m -p 0.025 -p 0.975 -o %s %s' % (out, file), shell=True)
