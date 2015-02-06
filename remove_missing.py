import gzip
import glob
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--chr", help="chromosome for which to run analysis")
args = parser.parse_args()
chr = args.chr

file = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/all_vcf/for_shapeit/gatk.ug.finch21.%s.allfilters.recoded_biallelicSNPs.vcf.gz' % chr

f = gzip.open(file, 'r')
out = file.replace('.vcf.gz', '.nomissing.vcf.gz')
o = gzip.open(out, 'w')
for l in f:
	if re.search('^#', l):
		o.write(l)
	else:
		count = len(re.findall('\.\/', l)) + len(re.findall('\/\.', l))
		if count != 42:
			o.write(l)
f.close()
o.close()
