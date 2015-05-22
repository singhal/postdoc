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
if sp == 'DBF':
	vcf_file = '/mnt/gluster/home/sonal.singhal1/DBF/after_vqsr/DB.allchrs.vqsr.filtered.snps.indels.vcf.gz'

out = '/mnt/lustre/home/sonal.singhal1/data/%s.nofilters.csv' % sp

f = gzip.open(vcf_file, 'r')
var = {}
for l in f:
	if re.search('PASS', l):
		d = re.split('\t', l.rstrip())

		# define if snp or indel
		type = 'snp'
		alleles = re.split(',', d[4]) + [d[3]]

		for allele in alleles:
			if len(allele) > 1:
				type = 'indel'

		# get the genotypes
		genos = []
		for geno in d[9:]:
			geno = re.search('^([^:]+)', geno).group(1)
			genos += re.split('[|/]', geno)
			
		# define if fixed or polymorphic	
		num_alleles = 0
		for ix, allele in enumerate(alleles):
			if genos.count(str(ix)) > 0:
				num_alleles += 1
		fixed_poly = 'fixed'
		if num_alleles > 1:
			fixed_poly = 'polymorphic'
		
		if d[0] not in var:
			var[d[0]] = {'snp': {'fixed': 0, 'polymorphic': 0}, 'indel': {'fixed': 0, 'polymorphic': 0}}

		var[d[0]][type][fixed_poly] += 1
f.close()

o = open(out, 'w')
o.write('chr,snp_indel,fixed_poly,count\n')
for chr in var:
	for type1 in var[chr]:
		for type2 in var[chr][type1]:
			o.write('%s,%s,%s,%s\n' % (chr, type1, type2, var[chr][type1][type2]))
o.close()
