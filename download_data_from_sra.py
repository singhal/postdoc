import re
import glob
import subprocess

files = glob.glob('*fastq')

for ix, file in enumerate(files):
	out = 'sra%s.sh' % ix
	o = open(out, 'w')
	#o.write('~/bin/sratoolkit.2.3.5-2-ubuntu64/bin/fastq-dump -I --split-files /mnt/gluster/home/sonal.singhal1/ficedula/sra_files/%s\n' % file)

	#fastq1 = file.replace('.sra', '_1.fastq')
	#fastq2 = file.replace('.sra', '_2.fastq')

	o.write('gzip /mnt/gluster/home/sonal.singhal1/ficedula/sra_files/%s\n' % file)
	#o.write('gzip /mnt/gluster/home/sonal.singhal1/ficedula/sra_files/%s\n' % fastq2)
	o.close()
