import glob
import subprocess
import re

files = glob.glob('/mnt/gluster/home/sonal.singhal1/*bam')

for sam in files:
	out = sam.replace('.bam', '.sorted')
	name = re.search('(SRR\d+)', sam).group(1)
	print "~/bin/samtools-0.1.19/samtools sort %s -o %s" % (sam, out)
	# subprocess.call('echo \"~/bin/samtools-0.1.19/samtools sort %s -o %s -m 28G\" | qsub -l h_vmem=34g -cwd -V -j y -N %s' % (sam, out, name), shell=True)
