import re
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument('--sp', help='species for which to run the analysis')
args = parser.parse_args()
chr = 'chrZ'
sp = args.sp

if sp == 'ZF':
	vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/all_vcf/gatk.ug.unrel_zf.chrZ.coverage.repeatmasked.filtered.nomendel.shared.recodedsex.vqsr2.vcf.gz'
	haps = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/finch19/%s_haplotypes.haps' % (chr)
elif sp == 'LTF':
	vcf = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.chrZ.coverage.repeatmasked.filtered.recodedsex.vqsr2.vcf.gz'
	haps = '/mnt/gluster/home/sonal.singhal1/LTF/phasing/PIR_approach/%s_haplotypes.haps' % (chr)

##############

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

out = vcf.replace('vqsr2.vcf.gz', 'phased.vqsr2.vcf.gz')
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
			tracker = 0
			for ix, geno in enumerate(d[9:]):
				geno = re.search('^([^:]+)', geno).group(1)
				if re.search('\/', geno):
					phase1 = allele_num[hap[pos][tracker]]
					phase2 = allele_num[hap[pos][tracker + 1]]
					
					d[ix+9] = re.sub(geno, '%s|%s' % (phase1, phase2), d[ix+9]) 
					tracker = tracker + 2
			hap.pop(pos, None)
			o.write('\t'.join(d) + '\n')
		else:
			# note that some of these sites are all fixed, so could be phased
			o.write(l)
f.close()
o.close()
