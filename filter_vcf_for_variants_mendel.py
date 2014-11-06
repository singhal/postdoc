import re
import gzip
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--chr", help="chromosome for which to run analysis")
args = parser.parse_args()
chr = args.chr

vcffile = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.%s.coverage.repeatmasked.vqsr.vcf.gz' % chr
vcfout = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.%s.coverage.repeatmasked.filtered.nomendel.vcf.gz' % chr

errorfile = '/mnt/gluster/home/sonal.singhal1/ZF/mendelian_errors/plink_results/all_zf.me.%s.mendel' % chr
sites_file = '/mnt/gluster/home/sonal.singhal1/ZF/mendelian_errors/all_zf.%s.map' % chr

sites_f = open(sites_file, 'r')
sites = {}
for l in sites_f:
	d = re.split('\s+', l)
	sites[d[1]] = int(d[3])
sites_f.close()

errors = {}
error_f = open(errorfile, 'r')
junk = error_f.next()
for l in error_f:
	d = re.split('\s+', l)
	errors[sites[d[4]]] = 1
error_f.close()

f = gzip.open(vcffile, 'r')
o = gzip.open(vcfout, 'w')

for l in f:
	if re.match('#', l):
		o.write(l)
	else:
		if re.search('PASS', l):
			d = re.split('\t', l)
			if int(d[1]) not in errors:
				o.write(l)
o.close()
f.close()
