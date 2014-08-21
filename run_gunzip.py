import glob
import subprocess

files = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/*gz')

for ix, file in enumerate(files):
	subprocess.call('echo "gunzip %s" | qsub -l h_vmem=1g -cwd -V -j y -N "gunz_%s"' % (file, ix), shell=True)
