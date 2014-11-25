import argparse
import glob
import gzip
import re

parser = argparse.ArgumentParser()
parser.add_argument("--sp", help="species for which to run analysis")
args = parser.parse_args()
sp = args.sp

if sp == 'ZF':
	vcfs = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/*phased*')
if sp == 'LTF':
	vcfs = glob.glob('/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/*phased*')
out = '/mnt/gluster/home/sonal.singhal1/%s/analysis/multiallelic_sites.csv' % sp
o = open(out, 'w')
o.write('chr,total_snps,multiallelic_snps\n')

for vcf in vcfs:
	all_snps = 0
	multiallele = 0

	chr = re.search('(chr[A-Z|0-9|a-z]+)', vcf).group(1)
	f = gzip.open(vcf, 'r')
	for l in f:
		if not re.search('#', l):
			d = re.split('\t', l.rstrip())
			alleles = [d[3]] + re.split(',', d[4])			
			ac = {}
			for ix, allele in enumerate(alleles):
				ac[str(ix)] = 0

			indels = False
			for allele in alleles:
				if len(allele) > 1:
					indels = True
		
			if not indels:
				for geno in d[9:]:
					geno = re.match('^([^:]+)', geno).group(1)
					geno = re.split('[|/]', geno)
					for g in geno:
						if g != '.':
							ac[g] += 1
				
				num_alleles = 0
				for allele, count in ac.items():
					if count > 0:
						num_alleles += 1
				
				if num_alleles > 1:
					all_snps += 1
					if num_alleles > 2:
						multiallele += 1
	f.close()
	o.write('%s,%s,%s\n' % (chr, all_snps, multiallele))
o.close()
					
