import re
import glob
import sys
import gzip

# directory with the VCF files that need to be cleaned up
dir = '/mnt/lustre/home/sonal.singhal1/LTF/after_vqsr/by_chr/'
# glob out those VCFs that still have "bad" sites
vcfs = []
for vcf in glob.glob(dir + '*vcf.gz'):
	if not re.search("biallelicSNPs", vcf):
		vcfs.append(vcf)

for vcf in vcfs:
	out = vcf
	out = out.replace('coverage', 'coverage.biallelicSNPs')
	i = gzip.open(vcf, 'r')
	o = gzip.open(out, 'w')

	for l in i:
		if re.search('^#', l):
			o.write(l)
		else:
			d = re.split('\t', l)
			if len(d[3]) == 1 and len(d[4]) == 1:
				o.write(l)
	i.close()
	o.close()

