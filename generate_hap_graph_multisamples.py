import random
import subprocess
import glob
import re

files = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/phase_samples_finch19/*graph')

for file in files:
	chr = re.search('(chr[A-Z|0-9]+)', file).group(1)
	for i in range(3):
		seed = random.randint(0,1000)
		out = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/phase_samples_finch19/%s_haplotypes.%s' % (chr, i)
		subprocess.call('/mnt/lustre/home/sonal.singhal1/bin/shapeit_v2r790/shapeit -convert --input-graph %s --output-sample %s --seed %s' % (file, out, seed), shell=True)
	
