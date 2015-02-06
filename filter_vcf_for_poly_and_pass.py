import re
import glob
import gzip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--chr", help="chromosome for which to run analysis")
args = parser.parse_args()
chr = args.chr

vcf = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.%s.coverage.repeatmasked.vqsr2.vcf.gz' % (chr)
out = vcf.replace('repeatmasked', 'repeatmasked.filtered')

out_f = gzip.open(out, 'w')
in_f = gzip.open(vcf, 'r')
	
for l in in_f:
	if re.search('^#', l):
		out_f.write(l)
	else:
		d = re.split('\t', l)
		if d[4] != '.':
			if d[6] == 'PASS':
				out_f.write(l)
out_f.close()
in_f.close()
