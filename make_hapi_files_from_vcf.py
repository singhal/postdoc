import re
import glob
import argparse
import gzip
import sys
import numpy

parser = argparse.ArgumentParser(description='Make HAPI files from VCF.')
parser.add_argument('--chr', help='chromosome')
args = parser.parse_args()
chr = args.chr

dir = '/mnt/gluster/home/sonal.singhal1/ZF/'
vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/gatk.ug.all_zf.%s.coverage.filtered.repeatmasked.vqsr.vcf.gz' % chr
geno_out_file = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/%s.hapi.geno' % chr
sites_out_file = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/%s.hapi.sites' % chr
marker_out_file = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/%s.hapi.list' % chr

ind_info = {	'MP1': ['fam1', '10', 'x', 'x', '2', '0'],
		'MP2': ['fam1', '12', 'x', 'x', '1', '0'],
		'MP3': ['fam1', '21', '12', '10', '1', '0'],
		'MP4': ['fam1', '22', '12', '10', '1', '0'],
		'MP5': ['fam1', '23', '12', '10', '1', '0']}

genos = {}
inds = []

errorfile = dir + 'mendelian_errors/plink_results/all_zf.me.%s.mendel' % chr
sites_file = dir + 'mendelian_errors/all_zf.%s.map' % chr

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
sites = {}

v = gzip.open(vcf, 'r')
for l in v:
	if re.search('^#CHROM', l):
		d = re.split('\t', l.rstrip())
		inds = d[-5:]
	if not re.search('^#', l):
		d = re.split('\t', l.rstrip())
		alleles = [d[3]] + re.split(',', d[4])
		indel = False
		for allele in alleles:
			if len(allele) > 1:
				indel = True
		if not indel and len(alleles) == 2:
			if int(d[1]) not in errors:
				all_missing = True
				for geno in d[-5:]:
					geno = re.search('^(\S/\S)', geno).group(1)
					if geno != './.':
						all_missing = False
				if not all_missing:
					genos[int(d[1])] = {}
					for ind, geno in zip(inds, d[-5:]):
						geno1 = re.search('^(\S)', geno).group(1)
						geno2 = re.search('\/(\S)', geno).group(1)
						if geno1 == '.':
							geno1 = 0
						else:
							geno1 = int(geno1) + 1
						if geno2 == '.':
							geno2 = 0
						else:
							geno2 = int(geno2) + 1
						genos[int(d[1])][ind] = '%s/%s' % (geno1, geno2)
v.close()

numbers = {'chr23': 1, 'chr24': 1, 'chr25': 1, 'chr26': 1, 'chr27': 1, 'chr28': 1, 'chr1A': 1, 
		'chr1B': 1, 'chr4A': 1, 'chrLG2': 1, 'chrLG5': 1, 'chrLGE22': 1, 'chrZ': 1}
sites_out_f = open(sites_out_file, 'w')
marker_out_f = open(marker_out_file, 'w')
chr_no = None
if chr in numbers:
	chr_no = numbers[chr]
else:
	chr_no = re.search('(\d+)', chr).group(1)
max_site = numpy.max(genos.keys())
for ix, site in enumerate(sorted(genos.keys())):
	marker_out_f.write('M rs%s\n' % ix)
	cm = site * 50 / float(max_site)
	sites_out_f.write('%s rs%s %.6f\n' % (chr_no, ix, cm))
sites_out_f.close()
marker_out_f.close()

geno_out_f = open(geno_out_file, 'w')
for ind in inds:
	geno_out_f.write('\t'.join(ind_info[ind]))
	for site in sorted(genos.keys()):
		geno_out_f.write('\t' + genos[site][ind])
	geno_out_f.write('\n')
geno_out_f.close()
