import glob
import subprocess

#' ~/bin/samtools-0.1.19/

dir = '/mnt/gluster/home/sonal.singhal1/LTF/bams/'
bams = glob.glob(dir + '*bam')

for ix, bam in enumerate(bams):
	print '%s\t%s' % (ix, bam)
	call = '~/bin/samtools-0.1.19/samtools index %s' % bam
	subprocess.call('echo \"%s\" | qsub -cwd -V -j y -N \'index%s\' -l h_vmem=8g' % (call, ix), shell=True)

