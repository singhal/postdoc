import glob
import subprocess

files = glob.glob('/mnt/lustre/home/sonal.singhal1/Darwin/g_fortis/SRR*sra')

for file in files:
	subprocess.call('~/bin/sratoolkit.2.3.5-2-ubuntu64/bin/fastq-dump %s -O /mnt/lustre/home/sonal.singhal1/Darwin/g_fortis/ --gzip' % file, shell=True)
	subprocess.call('rm %s' % file, shell=True)
