import re
import gzip
import argparse
import glob

parser = argparse.ArgumentParser()
parser.add_argument("--sp", help="species for which to run analysis")
args = parser.parse_args()
sp = args.sp

if sp == 'ZF':
	files = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/*coverage.repeatmasked.vqsr2*')
if sp == 'LTF':
	files = glob.glob('/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/*phased*')
out = '/mnt/lustre/home/sonal.singhal1/data/%s.coverage_repeat.csv' % sp

var = {}
for vcf_file in files:
	f = gzip.open(vcf_file, 'r')
	for l in f:
		if re.search('PASS', l):
			d = re.split('\t', l.rstrip())
			type = 'snp'
			alleles = re.split(',', d[4]) + [d[3]]
			for allele in alleles:
				if len(allele) > 1:
					type = 'indel'
			if d[0] not in var:
				var[d[0]] = {'snp': 0, 'indel': 0}
			var[d[0]][type] += 1
f.close()

o = open(out, 'w')
o.write('chr,snp,indel\n')
for chr in var:
	o.write('%s,%s,%s\n' % (chr, var[chr]['snp'], var[chr]['indel']))
o.close()
