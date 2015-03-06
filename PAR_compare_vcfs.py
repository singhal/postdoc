import re
import gzip

file1 = '/mnt/gluster/home/sonal.singhal1/PAR/zf/zf.cortex.filtered.chrZ_random.vcf'
file2 = '/mnt/gluster/home/sonal.singhal1/PAR/zf/gatk.ug.unrel_zf.chrZ_random.raw.snps.vcf'
out = '/mnt/gluster/home/sonal.singhal1/PAR/zf/unrel_zf.gatk_cortex_insersection.trusted_snps.vcf'

snps = {}
f = open(file1, 'r')
for l in f:
	if not re.search('^#', l):
		d = re.split('\t', l.rstrip())
		
		indel = False
		alleles = [d[3]] + re.split(',', d[4])
		for allele in alleles:
			if len(allele) > 1:
				indel = True

		if not indel:
			snps[d[1]] = alleles
f.close()

o = open(out, 'w')
f = open(file2, 'r')
for l in f:
	if re.search('^#', l):
		o.write(l)
	else:
		d = re.split('\t', l.rstrip())
		if d[1] in snps:
			alleles2 = [d[3]] + re.split(',', d[4])
			alleles1 = snps[d[1]]
			
			match = 0
			for allele in alleles2:
				if allele in alleles1:	
					match += 1
		
			if match >= 2:
				o.write(l)
f.close()

