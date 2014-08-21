import glob
import sys
import re
import gzip

min_allele = 1
dir = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/'
out = '/mnt/gluster/home/sonal.singhal1/ZF/segregating_sites_af%s.txt' % min_allele
vcfs = glob.glob(dir + '*biallelic*.vcf.gz')

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
				allele_counts = dict()
				for i, allele in enumerate(alleles):
					count = len(re.findall(str(i) + "/", line)) + len(re.findall("/" + str(i), line))
					if count > min_allele:
						allele_counts[i] = count		

				if len(allele_counts) == 2:
					seg_sites += 1
	o.write("%s\t%s\n" % (vcf, seg_sites))
	infile.close()
o.close() 
					
