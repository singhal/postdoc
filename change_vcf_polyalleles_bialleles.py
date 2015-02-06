import re
import gzip
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--chr", help="chromosome for which to run analysis")
args = parser.parse_args()
chr = args.chr

vcf = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.%s.coverage.repeatmasked.filtered.vqsr2.vcf.gz' % chr
vcfout = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/for_shapeit/gatk.ug.ltf.%s.allfilters.recoded_biallelicSNPs.vcf.gz' % chr

# designed to chop off some individuals, so look at this carefully

infile = gzip.open(vcf, 'r')
outfile = gzip.open(vcfout, 'w')

for l in infile:
	if re.search('^#', l):
		if re.search('^#CHROM', l):
			d = re.split('\t', l.rstrip())
			outfile.write('\t'.join(d[:]) + '\n')
		else:
			outfile.write(l)
	else:
		d = re.split('\t', l.rstrip())

		indel = False
		alleles = [d[3]] + re.split(",", d[4])
		for allele in alleles:
			if len(allele) > 1:
				indel = True

		if not indel:
			actual_alleles = []
			for ix, allele in enumerate(alleles):
				count = len(re.findall('%s\/' % ix, l)) + len(re.findall('\/%s' % ix, l))
				if count > 0:
					actual_alleles.append(allele)

			if len(actual_alleles) == 2:
				if len(d[4]) == 1:
					outfile.write('\t'.join(d[:]) + '\n')
				else:
					tmp = l.rstrip()
					for ix, allele in enumerate(alleles):
						if allele in actual_alleles:
							recode = actual_alleles.index(allele)
							tmp = tmp.replace('%s/' % ix, '%s/' % recode)
							tmp = tmp.replace('/%s' % ix, '/%s' % recode)	
					d2 = re.split('\t', tmp)

					outfile.write('\t'.join(d2[0:3]) + '\t' + actual_alleles[0] + 
							'\t' + actual_alleles[1] + '\t.\t' + d2[6] 
							+ '\t.\t' +  '\t'.join(d2[8:]) + '\n')
infile.close()
outfile.close()
				
