import operator
import re
import glob
import os.path
import gzip

# the minimum count at which sites will be dropped; i.e., if min_freq = 1, all singletons are dropped
min_freq = 1
files = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/results/*LGE22*_haplotypes.haps')
out_file = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/comparephasing_minfreq%s.txt' % min_freq
out = open(out_file, 'w')

def compare_ab(a,b):
	if a != '-' and b != '-':
		if a != b:
			return 'F'
		elif a == b:
			return 'T'
	else:
		return '0'


def compare_haplo(hap1, hap2):
	compare = map(compare_ab, hap1, hap2)
	compare = filter(lambda x: x != '0', compare)
                  
	bk_pts = 0
	for x, y in zip(compare, compare[1:]):
		if x != y:
			bk_pts += 1
	return bk_pts


for file1 in files:
	chr = re.search('chr([a-z|A-Z|0-9]+)', file1).group(1)
	file2 = file1.replace('PIR', 'family')

	site_freq = {}
	vcf_file = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/gatk.ug.zf.chr%s.filtered.coverage.recoded_biallelic.vqsr.vcf.gz' % chr
	vcf = gzip.open(vcf_file, 'r')
	for l in vcf:
		if not re.search('^#', l):
			d = re.split('\t', l)
			allele1 = len(re.findall('0\/', l)) + len(re.findall('\/0', l))
			allele2 = len(re.findall('1\/', l)) + len(re.findall('\/1', l))
			min_allele = min(allele1, allele2)
			
			site_freq[int(d[1])] = min_allele
	vcf.close()

	pir = {}
	site_index = {}
	pir_file =  '/mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/results/chr%s_PIRlist' % chr
	pirf = open(pir_file, 'r')
	num_sites = int(re.search('MAP\s+(\d+)', pirf.next()).group(1))
	for i in range(num_sites):
		site = re.search('^(\d+)', pirf.next()).group(1)
		site_index[i] = site
	for l in pirf:
		for match in re.findall('(\d+)\s+[A|T|C|G]', l):
			pir[site_index[int(match)]] = 1
	pirf.close()	

	if os.path.isfile(file2):	
		# these are the columns in the haplotype file that represent the
		# 	individuals we want to sample
		# Do not want to sample family zf
		# Only sampling one haplotype from each individual
		inds = range(5,43,2)

		haps1_a = {}
		haps1_b = {}
		sites = {}

		# initialize the dictionary
		for ind in inds:
			haps1_a[ind] = ''
			haps1_b[ind] = ''
			sites[ind] = list()

		f1 = open(file1, 'r')
		for l in f1:
			d = re.split('\s+', l)
			
			# create the haplotype
			for ind in inds:
				if d[ind] != d[ind+1]:
					if site_freq[int(d[2])] > min_freq:	
						if d[2] in pir:	
							sites[ind].append( int(d[2]) )
							haps1_a[ind] += d[ind]
							haps1_b[ind] += d[ind + 1]
		haps2_tmp = {}
		haps2 = {}
		for ind in inds:
			haps2_tmp[ind] = dict()
			haps2[ind] = ''

		f2 = open(file2, 'r')
		for l in f2:
			d = re.split('\s+', l)
			
			# create the haplotype
			for ind in inds:
				haps2_tmp[ind][ int(d[2]) ] = d[ind]

		for ind in inds:
			for site in sites[ind]:
				if site in haps2_tmp[ind]:
					haps2[ind] += haps2_tmp[ind][site]
				else:
					haps2[ind] += '-'

		for ind in inds:
			bk_pts1 = compare_haplo(haps1_a[ind], haps2[ind])
			bk_pts2 = compare_haplo(haps1_b[ind], haps2[ind])

			bk_pts = min(bk_pts1, bk_pts2)

			diff = bk_pts / float( len(sites[ind]) )

			out.write('%s\t%s\t%s\t%.3f\n' % (chr, ind, len(sites[ind]), diff))
