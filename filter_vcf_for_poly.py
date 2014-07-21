import re
import glob
import gzip

dir1 = '/mnt/lustre/home/sonal.singhal1/ZF/phasing/family_approach/family_vcf/by_chr/'
append1 = dir1 + 'gatk.ug.zf_fam.'
dir2 = '/mnt/lustre/home/sonal.singhal1/ZF/phasing/family_approach/unrel_vcf/by_chr/'
append2 = dir2 + 'gatk.ug.zf_unrel.'

chrs = ['chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr1A', 'chr1B', 
	'chr1', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28',
	'chr2', 'chr3', 'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chrLG2', 'chrLG5', 'chrLGE22', 'chrZ']

def variants_in_vcf(vcffile, var):
	v = gzip.open(vcffile, 'r')
	for l in v:
		if re.search('PASS', l):
			d = re.split('\t', l)
			var[d[1]] = 1
	v.close()
	return var

def filter_sites_in_vcf(vcffile, var, filestart, chr):
	out = filestart + chr + '.coverage.no_mendel.poly_in_any_zf.vqsr.vcf.gz'	

	v = gzip.open(vcffile, 'r')
	o = gzip.open(out, 'w')
	
	for l in v:
		if re.search('^#', l):
			o.write(l)
		else:
			d = re.split('\t', l)
			if d[1] in var:
				o.write(l)

	o.close()
	v.close()
	
for chr in chrs:
	vcf1 = append1 + chr + '.coverage.no_mendel.vqsr.vcf.gz'
	vcf2 = append2 + chr + '.coverage.no_mendel.vqsr.vcf.gz'

	# identify variants in one or the other VCF
	var = dict()
	var = variants_in_vcf(vcf1, var)
	var = variants_in_vcf(vcf2, var)

	filter_sites_in_vcf(vcf1, var, append1, chr)
	filter_sites_in_vcf(vcf2, var, append2, chr)

	
