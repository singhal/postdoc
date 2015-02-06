import re
import gzip

# ped = '/mnt/gluster/home/sonal.singhal1/ZF/zf_inds_plink.ped'
ped = '/mnt/gluster/home/sonal.singhal1/LTF/LTF_inds_plink.ped'
# vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/all_vcf/for_shapeit/gatk.ug.finch19.chrZ.allfilters.recodedsex.recoded_biallelicSNPs.vcf.gz'
vcf = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/for_shapeit/gatk.ug.ltf.chrZ.allfilters.recodedsex.recoded_biallelicSNPs.vcf.gz'

# output data
out = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/for_shapeit/gatk.ug.ltf.chrZ.allfilters.recodedsex.recoded_biallelicSNPs.males.vcf.gz'
# out = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/all_vcf/for_shapeit/gatk.ug.finch19.chrZ.allfilters.recodedsex.recoded_biallelicSNPs.males.vcf.gz'
# outdir = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/all_vcf/for_shapeit/'
outdir = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/for_shapeit/'

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
snps = {}
num_fem = 0

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
		
		newline = '\t'.join(d[0:9])
	        males = []
		fem_gens = []
                genos = []
		for ind, geno in zip(inds, d[9:]):
                       	if sex[ind] == 'M':
                       	        males.append(geno)
				geno = re.search('(\S\/\S)', geno).group(1)
				genos += re.split('/', geno)
			else:
				allele = re.search('^(\S)', geno).group(1)
				fem_gens.append(allele)

		if len(genos) > genos.count('.'):
			num_fem = len(fem_gens)
			double = []
			for gen in fem_gens:
				if gen == '.':
					double.append('?')
					double.append('?')
				else:
					double.append(gen)
					double.append(gen)
			snps[int(d[1])] = {'ref': d[3], 'alt': d[4], 'gen': double}

        	        newline = newline + '\t' + '\t'.join(males)
        	        o.write(newline + '\n')
o.close()
f.close()

samplefile = '%schrZ.sample' % outdir
o = open(samplefile, 'w')
o.write('sample population group sex\n')
for ix in range(num_fem):
        o.write('female%s female%s female%s 2\n' % (ix, ix, ix))
o.close()

legendfile = '%schrZ.legend.gz' % outdir
o = gzip.open(legendfile, 'w')
o.write('id position a0 a1\n')
for ix, pos in enumerate(sorted(snps.keys())):
        o.write('SNP%s %s %s %s\n' % (ix, pos, snps[pos]['ref'], snps[pos]['alt']))
o.close()

hapfile = '%schrZ.hap.gz' % outdir
o = gzip.open(hapfile, 'w')
for pos in sorted(snps.keys()):
        o.write(' '.join(snps[pos]['gen']) + '\n')
o.close()
					
