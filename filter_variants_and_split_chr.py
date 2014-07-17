import gzip
import re

vcf_file = '/mnt/gluster/home/emleffler/genotype_callsets/zebrafinch/zf_family/unified_genotyper/after_vqsr/gatk.ug.MP1-5.allchrs.snps.indels.vqsr2.allsites.vcf.gz'

outdir = '/mnt/lustre/home/sonal.singhal1/ZF/phasing/family_approach/family_vcf/'

chr = 'NA'
mainfile = gzip.open(vcf_file, 'r')
header = []
for l in mainfile:
	if re.search('^#', l):
		header.append(l)
	else:
		d = re.split('\t', l)
		if d[0] != chr:
			if chr != 'NA':
				out.close()
			chr = d[0]
			outfile = outdir + 'gatk.ug.MP1-5.%s.snps.indels.variable.vqsr.vcf.gz' % chr
			out = gzip.open(outfile, 'w')
			for l in header:
				out.write(l)
		if re.search('AF=\d', l):
			out.write(l)
mainfile.close()

