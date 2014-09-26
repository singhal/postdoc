import glob
import re
import subprocess
import os

chrs = ['chr1', 'chr1A', 'chr1B', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
	'chr10', 'chr11', 'chr12', 'chr4A', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
	'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28' ]

#chrs = ['chr16', 'chrLGE22']

haplo1 = 'haplo0'
haplo2 = 'haplo39'
out1 = '/mnt/gluster/home/sonal.singhal1/LTF/phasing/%s_tmp.fa' % haplo1
out2 = '/mnt/gluster/home/sonal.singhal1/LTF/phasing/%s_tmp.fa' % haplo2
outfile = '/mnt/gluster/home/sonal.singhal1/LTF/analysis/LTF.coalhmm.fa'

haplo1seq = ''
haplo2seq = ''

for chr in chrs:
	seq = '/mnt/gluster/home/sonal.singhal1/LTF/phasing/PIR_approach/%s_haplotypes.fasta' % chr
	subprocess.call("~/bin/samtools-0.1.19/samtools faidx %s" % seq, shell=True)
	subprocess.call("~/bin/samtools-0.1.19/samtools faidx %s %s > %s" % (seq, haplo1, out1), shell=True)
	subprocess.call("~/bin/samtools-0.1.19/samtools faidx %s %s > %s" % (seq, haplo2, out2), shell=True)

	f1 = open(out1, 'r')
	f2 = open(out2, 'r')

	for l1, l2 in zip(f1, f1):
		if not re.search('>', l1):
			l1 = list(l1.rstrip())
			l2 = list(l2.rstrip())

			for bp1, bp2 in zip(l1, l2):
				if bp1 != 'N' and bp2 != 'N':
					haplo1seq += bp1
					haplo2seq += bp2

	f1.close()
	f2.close()

out = open(outfile, 'w')
out.write('>%s\n%s\n' % (haplo1, haplo1seq))
out.write('>%s\n%s\n' % (haplo2, haplo2seq))			
out.close()

os.remove(out1)
os.remove(out2)
