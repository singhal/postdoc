import glob
import subprocess
import re

files = glob.glob("*.fa")

for ix, file in enumerate(files):
	theta = re.search('theta([0-9|\.]+)', file).group(1)
	subprocess.call("/mnt/lustre/home/sonal.singhal1/bin/LDhelmet_v1.6/ldhelmet find_confs -w 50 -o out_theta%s_%s %s" % (theta, ix, file), shell=True)  
