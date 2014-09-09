import glob
import subprocess

files = glob.glob('/mnt/gluster/home/sonal.singhal1/*F/analysis/PSMC/*psmcfa')
for ix, file in enumerate(files):
	out = file.replace('.psmcfa', '.psmc')
	call = '~/bin/psmc/psmc -N25 -t10 -r5 -p \"4+25*2+4+6\" -o %s %s' % (out, file)
	subprocess.call('echo \'%s\' | qsub -l h_vmem=20g -cwd -V -j y -N \'psmc%s\'' % (call, ix), shell=True)

