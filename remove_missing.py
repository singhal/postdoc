import gzip
import glob
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--chr", help="chromosome for which to run analysis")
args = parser.parse_args()
chr = args.chr

file = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/gatk.ug.all_zf.%s.coverage.filtered.repeatmasked.recoded_biallelicSNPs.nomendel.vcf.gz' % chr

f = gzip.open(file, 'r')
out = file.replace('nomendel', 'nomendel.nomissing')
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
