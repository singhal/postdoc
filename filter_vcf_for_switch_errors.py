import re
import glob
import gzip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--chr", help="chromosome for which to run analysis")
args = parser.parse_args()
chr = args.chr

vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/all_vcf/gatk.ug.all_zf.%s.coverage.repeatmasked.filtered.nomendel.shared.vqsr2.vcf.gz' % (chr)
out = vcf.replace('shared', 'shared.noswitch')
switches = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/breaks/badintervals.%s.hapi.csv' % chr

ditch = {}
f = open(switches, 'r')
header = f.next()
for l in f:
	d = re.split(',', l.rstrip())
	for i in range(int(d[0])+1, int(d[1])):
		ditch[i] = 1
f.close()

out_f = gzip.open(out, 'w')
in_f = gzip.open(vcf, 'r')
	
for l in in_f:
	if re.search('^#', l):
		out_f.write(l)
	else:
		d = re.split('\t', l)
		if int(d[1]) not in ditch:
			out_f.write(l)
out_f.close()
in_f.close()
