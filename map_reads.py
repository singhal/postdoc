import glob
import re
import subprocess

reads = glob.glob('/mnt/lustre/home/sonal.singhal1/Darwin/*gz')
dir = '/mnt/lustre/home/sonal.singhal1/Darwin/'

for ix, new_read in enumerate(reads):
	# new_read = dir + 'finch%s.fastq' % ix
	# subprocess.call('mv %s %s' % (read, new_read), shell=True)
	# subprocess.call('gzip %s' % (new_read), shell=True)
	# new_read += '.gz'
	job = dir + 'map%s.sh' % (ix)
	out = open(job, 'w')
	out.write('/mnt/lustre/home/sonal.singhal1/bin/stampy-1.0.23/stampy.py -g ~/reference/taeGut1.bamorder -h ~/reference/taeGut1.bamorder -M %s --substitutionrate=0.12 -o ~/Darwin/finch%s.sam --sanger' % (new_read, ix))
	out.close()
	subprocess.call('chmod a+x %s' % job, shell=True)
	subprocess.call('echo \"%s\" | qsub -l h_vmem=5g -cwd -V -j y -N \'finch%s\'' % (job, ix), shell=True)

