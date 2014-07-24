import glob
import subprocess
import re

# files = glob.glob('/mnt/lustre/home/sonal.singhal1/Darwin/*gz')
redo = ['10', '11', '4', '5', '7', '8']

for num in redo:
	name = 'finch' + str(num)
	file = '/mnt/lustre/home/sonal.singhal1/Darwin/' + name + '.fastq.gz'
	
	call = 'echo "~/bin/bowtie2-2.2.3/bowtie2 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 --local -p 8 -x ~/reference/taeGut1.bamorder.fasta -U %s -S %s.sam" | qsub -l h_vmem=6g -cwd -j y -V -N \'%s\'' % (file, name, name)

	subprocess.call('%s' % call, shell=True)
