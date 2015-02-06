import re
import glob
import gzip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--chr", help="chromosome for which to run analysis")
args = parser.parse_args()
chr = args.chr

vcf1 = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.%s.coverage.repeatmasked.filtered.nomendel.vqsr2.vcf.gz' % (chr)
vcf2 = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/all_vcf/gatk.ug.all_zf.%s.coverage.repeatmasked.filtered.nomendel.vqsr2.vcf.gz' % chr
out = vcf2.replace('nomendel', 'nomendel.shared')

sites = {}
f = gzip.open(vcf1, 'r')
for l in f:
        if not re.search('^#', l):
                d = re.split('\t', l)
                sites[ int(d[1]) ] = 1
f.close()

out_f = gzip.open(out, 'w')
in_f = gzip.open(vcf2, 'r')
	
for l in in_f:
	if re.search('^#', l):
		out_f.write(l)
	else:
		d = re.split('\t', l)
		if int(d[1]) in sites:
                        out_f.write(l)
out_f.close()
in_f.close()
