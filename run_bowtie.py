import glob
import subprocess
import re

files = glob.glob('/mnt/lustre/home/sonal.singhal1/Darwin/g_fortis/*gz')

for file in files:
	name = re.search('(SRR.*)\.fast', file).group(1)
	out = '/mnt/gluster/home/sonal.singhal1/' + name + '.sam'
	
	call = 'echo "~/bin/bowtie2-2.2.3/bowtie2 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 --local -p 8 -x ~/reference/taeGut1.bamorder.fasta -U %s -S %s" | qsub -l h_vmem=6g -cwd -j y -V -N \'%s\'' % (file, out, name)

	subprocess.call('%s' % call, shell=True)
