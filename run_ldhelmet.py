import glob
import re
import subprocess

files = glob.glob('/mnt/gluster/home/sonal.singhal1/LTF/phasing/PIR_approach/*fasta')

for file in files:
	chr = re.search('(chr\S+)_hap', file).group(1)
	out = '/mnt/gluster/home/sonal.singhal1/LTF/analysis/LDhelmet/%s.conf' % chr
	subprocess.call('echo \"~/bin/LDhelmet_v1.6/ldhelmet find_confs -w 50 --num_threads 12 -o %s %s\" | qsub -l h_vmem=10g -cwd -V -j y -N %s' % (out, file, chr), shell=True)
