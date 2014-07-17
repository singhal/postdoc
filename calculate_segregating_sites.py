import glob
import sys
import re
import gzip

dir = '/mnt/lustre/home/sonal.singhal1/ZF/after_vqsr/by_chr/'
out = '/mnt/lustre/home/sonal.singhal1/ZF/segregating_sites.txt'
vcfs = []
for vcf in glob.glob(dir + '*.vcf.gz'):
	if not re.search('biallelic', vcf):
		vcfs.append(vcf)

o = open(out, 'w')
for vcf in vcfs:
	seg_sites = 0
	infile = gzip.open(vcf, 'r')
	for line in infile:
		if not re.search('^#', line):
			d = re.split('\t', line)
			alleles = re.split(',', d[4])
			indel = False
			if len(d[3]) > 1:
				indel = True
			for allele in alleles:
				if len(allele) > 1:
					indel = True
			if not indel:
				num_alleles = 1 + len(alleles)
				allele_counts = dict()
				for i in range(num_alleles):
					allele_counts[i] = len(re.findall(str(i) + "/", line)) + len(re.findall("/" + str(i), line))
				num_poly = 0
				for allele in allele_counts:
					if allele_counts[allele] > 0:
						num_poly += 1
				if num_poly > 1:
					seg_sites += 1
	o.write("%s\t%s\n" % (vcf, seg_sites))
	infile.close()
o.close() 
					
