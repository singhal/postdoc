import re
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument('--sp', help='species for which to run the analysis')
parser.add_argument('--chr', help='chromosome for which to run the analysis')
args = parser.parse_args()
chr = args.chr
sp = args.sp

if sp == 'ZF':
	vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.%s.coverage.repeatmasked.filtered.nomendel.vcf.gz' % chr
elif sp == 'LTF':
	vcf = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.%s.filtered.coverage.repeatmasked.vqsr.vcf.gz' % chr

##############
haps = '/mnt/gluster/home/sonal.singhal1/%s/phasing/PIR_approach/results/%s_haplotypes.haps' % (sp, chr)

def snps(snp, ref, alt):
	if snp == '0':
		return ref
	else:
		return alt

hap = {}
f = open(haps, 'r')
for ix, l in enumerate(f):
	d = re.split('\s+', l.strip())
	pos = int(d[2])
	# translate 0 and 1 into the actual alleles
	hap[pos] = [snps(x, d[3], d[4]) for x in d[5:]]
f.close()

out = vcf.replace('vcf.gz', 'phased.vcf.gz')
f = gzip.open(vcf, 'r')
o = gzip.open(out, 'w')

for l in f:
	if re.match('#', l):
		o.write(l)
	else:
		d = re.split('\t', l.rstrip())
		pos = int(d[1])
		if pos in hap:
			# define the alleles by their number 
			alleles = [d[3]] + re.split(',', d[4])
			allele_num = {}
			for ix, allele in enumerate(alleles):
				allele_num[allele] = ix
			
			# figure out the genotype on each strand
			for ix, geno in enumerate(d[9:]):
				# should have probably written this as a tuple ..
				geno1 = allele_num[hap[pos][ix*2]]
				geno2 = allele_num[hap[pos][ix*2 + 1]]
				# replace the unphased genotype with the phased genotype
				d[ix+9] = re.sub('\S/\S', '%s|%s' % (geno1, geno2), d[ix+9])
			o.write('\t'.join(d) + '\n')
			# have to do this because some of the vcfs have two variants at the same position
			# one is an indel, the other is a SNP
			# SNPs always come first, and only want to take care of them.
			hap.pop(pos, None)
		else:
			# note that some of these sites are all fixed, so could be phased
			o.write(l)
