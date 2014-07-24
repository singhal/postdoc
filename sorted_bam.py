import glob
import subprocess
import re

files = glob.glob('/mnt/lustre/home/sonal.singhal1/Darwin/*sam')

for sam in files:
	out = sam.replace('.sam', '.sorted')
	name = re.search('(finch\d+)', sam).group(1)
	subprocess.call('echo \"samtools view -bS %s | samtools sort - %s\" | qsub -l h_vmem=5g -cwd -V -j y -N %s' % (sam, out, name), shell=True)
