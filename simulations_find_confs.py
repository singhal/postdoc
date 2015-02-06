import glob
import subprocess
import re

files = glob.glob("/mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/*.fasta")

for ix, file in enumerate(files):
	subprocess.call("/mnt/lustre/home/sonal.singhal1/bin/LDhelmet_v1.6/ldhelmet find_confs -w 50 -o out_%s %s" % (ix, file), shell=True)  
