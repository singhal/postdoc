import re
import gzip

ped = '/mnt/gluster/home/sonal.singhal1/ZF/zf_inds_plink.ped'
# ped = '/mnt/gluster/home/sonal.singhal1/LTF/LTF_inds_plink.ped'
vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.chrZ.coverage.repeatmasked.filtered.nomendel.vqsr2.vcf.gz'
# vcf = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.chrZ.coverage.repeatmasked.filtered.vqsr2.vcf.gz'
# output data
out = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.chrZ.coverage.repeatmasked.filtered.nomendel.recodedsex.vqsr2.vcf.gz'
# out = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.chrZ.coverage.repeatmasked.filtered.recodedsex.vqsr2.vcf.gz'
# proportion wrong; if 3 or more females are heterozygous at a given site, call it a bad site
filter = 3

#################################

# open out file
o = gzip.open(out, 'w')

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

inds = []
f = gzip.open(vcf, 'r')
for l in f:
	l = l.rstrip()
	if re.search('#CHROM', l):
		inds = re.split('\t', l)[9:]
		o.write(l + '\n')
	elif re.search('#', l):
		o.write(l + '\n')
	elif re.search('PASS', l):
		d = re.split('\t', l)
		miscalled_hets = 0
		for ind, geno in zip(inds, d[9:]):
			genomatch = re.search('(\S)\/(\S)', geno)
			het = False
			if genomatch.group(1) != genomatch.group(2):
				het = True
			# this is the error situation
			# females should never be hets because they only have one Z
			# (also, no PAR!)
			if sex[ind] == 'F' and het:
				miscalled_hets += 1
		if miscalled_hets < filter:
			for ix, (ind, geno) in enumerate(zip(inds, d[9:])):
				if sex[ind] == 'F':
					genomatch = re.search('(\S)\/(\S)', geno)
					geno1 = genomatch.group(1)
					geno2 = genomatch.group(2)
                        		het = False
                        		if geno1 != geno2:
                                		het = True
					if het:
						d[ix+9] = d[ix+9].replace('%s/%s' % (geno1, geno2), '.')
					else:
						d[ix+9] = d[ix+9].replace('%s/%s' % (geno1, geno2), geno1)
			o.write('\t'.join(d) + '\n')
