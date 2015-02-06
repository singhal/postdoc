import glob
import subprocess
import re

files = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/*vqsr2*')

for vcfgz in files:
	chr = re.search('(chr[A-Z|0-9|\_|a-z]+)', vcfgz).group(1)
	o = open('unrel_%s.sh' % chr, 'w')
	if re.search('gz', vcfgz):
		o.write('gunzip %s\n' % vcfgz)
		vcf = vcfgz.replace('.gz', '')
	else:
		vcf = vcfgz
	o.write('tail %s\n' % vcf)
	o.write('~/bin/tabix-0.2.6/bgzip %s\n' % vcf)
	o.write('~/bin/tabix-0.2.6/tabix -p vcf %s\n' % vcfgz)
	o.close()
