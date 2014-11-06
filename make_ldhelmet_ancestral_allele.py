import re
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--chr", help="chromosome for which to run")
parser.add_argument("--sp", help="species for which to run")
args = parser.parse_args()

chr = args.chr
sp = args.sp

nuc = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

for chr in [chr]:
	aa_file = '/mnt/gluster/home/sonal.singhal1/%s/ancestral_allele/ancestral_allele.%s.csv' % (sp, chr)
	sites_file = '/mnt/gluster/home/sonal.singhal1/%s/phasing/PIR_approach/%s_sites.csv' % (sp, chr)
	out_file = '/mnt/gluster/home/sonal.singhal1/%s/ancestral_allele/ancestral_allele.%s.ldhelmet.txt' % (sp, chr)

	anc_info = {}
	
	aa_f = open(aa_file, 'r')
	header = aa_f.next()
	for l in aa_f:
		d = re.split(',', l.rstrip())
		site = int(d[1])

		if d[-1] == 'N':
			anc_info[site] = [0.30, 0.20, 0.20, 0.30]
		else:
			anc_info[site] = [0.02, 0.02, 0.02, 0.02]
			anc_info[site][nuc[d[-1]]] = 0.94
	aa_f.close()

	out_f = open(out_file, 'w')
	sites_f = open(sites_file, 'r')
	for l in sites_f:
		d = re.split(',', l.rstrip())
		d[1] = int(d[1])
		if d[1] in anc_info:
			out_f.write('%s %s\n' % (d[1] - 1, ' '.join('%.2f' % i for i in anc_info[int(d[1])])))
		else:
			out_f.write('%s %s\n' % (d[1] - 1, '0.29 0.21 0.21 0.29'))
			print 'Undefined ancestral sequence: %s %s' % (chr, d[1])
	out_f.close()
	sites_f.close()
