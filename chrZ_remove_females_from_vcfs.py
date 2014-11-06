import re
import gzip

ped = '/mnt/gluster/home/sonal.singhal1/ZF/zf_inds_plink.ped'
# ped = '/mnt/gluster/home/sonal.singhal1/LTF/LTF_inds_plink.ped'
vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/for_shapeit/gatk.ug.unrel_zf.chrZ.coverage.repeatmasked.filtered.recodedsex.recoded_biallelicSNPs.nomendel.vcf.gz'
# vcf = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/for_shapeit/gatk.ug.ltf.chrZ.filtered.coverage.repeatmasked.recodedsex.recoded_biallelicSNPs.vcf.gz'
# output data
# out = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/for_shapeit/gatk.ug.ltf.chrZ.filtered.coverage.repeatmasked.recodedsex.recoded_biallelicSNPS.males.vcf.gz'
out = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/for_shapeit/gatk.ug.unrel_zf.chrZ.coverage.repeatmasked.filtered.recodedsex.recoded_biallelicSNPs.nomendel.males.vcf.gz'

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
		d = re.split('\t', l)
		inds = d[9:]
		newline = '\t'.join(d[0:9])
		males = []
		for ind in inds:
			if sex[ind] == 'M':
				males.append(ind)
		newline = newline + '\t' + '\t'.join(males)	
		o.write(newline + '\n')
	elif re.search('#', l):
		o.write(l + '\n')
	elif re.search('PASS', l):
		d = re.split('\t', l)
		genos = d[9:]
                newline = '\t'.join(d[0:9])
                males = []
                for ind, geno in zip(inds, genos):
                        if sex[ind] == 'M':
                                males.append(geno)
                newline = newline + '\t' + '\t'.join(males)
                o.write(newline + '\n')
o.close()
f.close()
					
