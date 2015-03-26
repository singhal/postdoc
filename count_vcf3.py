import re
import gzip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--sp", help="species for which to run analysis")
args = parser.parse_args()
sp = args.sp

if sp == 'ZF':
	vcf_file = '/mnt/gluster/home/emleffler/genotype_callsets/zebrafinch/zf_unrels/unified_genotyper/after_vqsr/gatk.ug.unrelzf.allchrs.snps.indels.vqsr2.vcf.gz'
if sp == 'LTF':
	vcf_file = '/mnt/gluster/home/emleffler/genotype_callsets/longtailedfinch/after_vqsr/gatk.ug.ltf.allchrs.snps.indels.vqsr2.vcf.gz'
out = '/mnt/lustre/home/sonal.singhal1/data/%s.nofilters.csv' % sp

f = gzip.open(vcf_file, 'r')
var = {}
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
