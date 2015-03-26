import glob
import subprocess
import re

files = glob.glob('/mnt/gluster/data/internal_restricted_supp/finches_2014/longtailedfinch/realigned_mateFixed_recal_bams/*bam')

files.append('/mnt/gluster/home/sonal.singhal1/DBF/bams/DB.allchrs.realigned.mateFixed.recal.bam')

for file in files:
	out = '/mnt/lustre/home/sonal.singhal1/data/' + re.sub('^.*bams\/', '', file) + '.flagstat'
	subprocess.call('samtools flagstat %s > %s' % (file, out), shell=True)

