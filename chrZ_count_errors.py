import re
import gzip

ped = '/mnt/gluster/home/sonal.singhal1/ZF/zf_inds_plink.ped'
# ped = '/mnt/gluster/home/sonal.singhal1/LTF/LTF_inds_plink.ped'
vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.chrZ.coverage.repeatmasked.vqsr.vcf.gz'
# vcf = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.chrZ.filtered.coverage.repeatmasked.vqsr.vcf.gz'
# output data
# proportion wrong; if 3 or more females are heterozygous at a given site, call it a bad site
filter = 3

#################################

# get sex data
sex = {}
f = open(ped, 'r')
for l in f:
	d = re.split('\s+', l)
	if d[4] == '1':
		sex[d[1]] = 'M'
	elif d[4] == '2':
		sex[d[1]] = 'F'
f.close()

all_female_geno = 0
wrong_female_geno = 0
filtered_sites = 0
all_sites = 0

inds = []
f = gzip.open(vcf, 'r')
for l in f:
	l = l.rstrip()
	if re.search('#CHROM', l):
		inds = re.split('\t', l)[9:]
	elif re.search('PASS', l):
		d = re.split('\t', l)
		all_sites += 1
		miscalled_hets = 0
		for ind, geno in zip(inds, d[9:]):
			genomatch = re.search('(\S)\/(\S)', geno)
			het = False
			if genomatch.group(1) != genomatch.group(2):
				het = True
			# this is the error situation
			# females should never be hets because they only have one Z
			# (also, no PAR!)
			if sex[ind] == 'F':
				if het:
					miscalled_hets += 1
					all_female_geno += 1
					wrong_female_geno += 1
				else:
					all_female_geno += 1
				
		if miscalled_hets >= filter:
			filtered_sites += 1

print 'all_sites: %s' % all_sites
print 'filtered_sites: %s' % filtered_sites
print 'all_female_geno: %s' % all_female_geno
print 'wrong_female_geno: %s' % wrong_female_geno
