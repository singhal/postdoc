import glob
import subprocess
import os
import re

files = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/maps/*bpen5')
files = filter(lambda x: not re.search('.txt', x), files)


for file in files:
	print file
	out = '%s.txt' % file
	if not os.path.isfile(out):
		subprocess.call('~/bin/LDhelmet_v1.6/ldhelmet post_to_text -p 0.025 -p 0.975 -m -o %s %s' % (out, file), shell=True)
