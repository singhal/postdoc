import glob
import sys
import re
import gzip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--sp", help="species [ZF|LTF]")
parser.add_argument("--min_allele", help="min allele count to be considered SS")
args = parser.parse_args()

sp = args.sp
min_allele = int(args.min_allele)

if sp == 'ZF':
	vcfs = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/*phased*')
if sp == 'LTF':
	vcfs = glob.glob('/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/*phased*')
out = '/mnt/gluster/home/sonal.singhal1/%s/analysis/segregatingsites_af%s.txt' % (sp, min_allele)

o = open(out, 'w')
for vcf in vcfs:
	seg_sites = 0
	infile = gzip.open(vcf, 'r')
	for line in infile:
		if not re.search('^#', line):
			d = re.split('\t', line)
			alleles = [d[3]] + re.split(',', d[4])
			indel = False
			for allele in alleles:
				if len(allele) > 1:
					indel = True
			if not indel:
				genos = []
				for geno in d[9:]:
					geno = re.search('^([^:]+)', geno).group(1)
                                        genos += re.split('[|/]', geno)
                                num_alleles = 0
                                for ix, allele in enumerate(alleles):
                                        if genos.count(str(ix)) > min_allele:
                                                num_alleles += 1
				if num_alleles > 1:
					seg_sites += 1
	o.write("%s\t%s\n" % (vcf, seg_sites))
	infile.close()
o.close()
