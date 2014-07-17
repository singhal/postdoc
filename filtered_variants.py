import gzip
import re

file = '/mnt/gluster/home/emleffler/genotype_callsets/zebrafinch/zf_family/unified_genotyper/after_vqsr/gatk.ug.MP1-5.allchrs.snps.indels.vqsr2.allsites.vcf.gz'
out = '/mnt/lustre/home/sonal.singhal1/ZF/phasing/family_approach/family_vcf/gatk.ug.MP1-5.allchrs.snps.indels.vqsr.filtered.vcf.gz'

r = gzip.open(file, 'r')
o = gzip.open(out, 'w')

for l in r:
	if re.search('^#', l):
		o.write(l)
	else:
		if re.search('AF=\d', l):
			if re.search('PASS', l):
				o.write(l)
o.close()
r.close()
