import glob
import re
import gzip
import random
from itertools import izip

file = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.%s.filtered.coverage.recoded_biallelicSNPs.vqsr.vcf.gz'
chrs = ['chr1', 'chr1A', 'chr1B', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr4A', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28'] 

out_file = '/mnt/gluster/home/sonal.singhal1/LTF/analysis/LTF.geno'
out = open(out_file, 'w')
subsample = 100

for chr in chrs:
	vcf = file % chr
	f = gzip.open(vcf)
	for ix, l in enumerate(f):
		if not ix % subsample:
			if not re.search('^#', l):
				snp = ''
				d = re.split('\t', l.rstrip())
				for geno in d[9:]:
					if re.match('0/0', geno):
						snp += '0'
					elif re.match('0/1', geno):
						snp += '1'
					elif re.match('1/1', geno):
						snp += '2'
					else:
						snp += '9'
				out.write(snp + '\n')
	f.close()
out.close()
