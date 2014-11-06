import re
import gzip

# vcffile = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/for_shapeit/gatk.ug.ltf.chrZ.filtered.coverage.repeatmasked.recodedsex.recoded_biallelicSNPs.vcf.gz'
# outdir = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/for_shapeit/'
# missing = [12297573, 25324470, 31785435, 61311164]

vcffile = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/for_shapeit/gatk.ug.unrel_zf.chrZ.coverage.repeatmasked.filtered.recodedsex.recoded_biallelicSNPs.nomendel.vcf.gz'
outdir = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/for_shapeit/'
missing = [15901058]

f = gzip.open(vcffile, 'r')

num_fem = 10
snps = {}

for l in f:
	if not re.search('#', l):
		l = l.rstrip()
		d = re.split('\t', l)
		fem_gens = []
		for geno in d[9:]:
			if not re.search('\S\/', geno):
				allele = re.search('^(\S)', geno).group(1)
				fem_gens.append(allele)
		
		if int(d[1]) not in missing:
			double = []
			for gen in fem_gens:
				if gen == '.':
					double.append('?')
					double.append('?')
				else:
					double.append(gen)
					double.append(gen)
				snps[int(d[1])] = {'ref': d[3], 'alt': d[4], 'gen': double}
f.close()

samplefile = '%sfemalehaps.sample' % outdir
o = open(samplefile, 'w')
o.write('sample population group sex\n')
for ix in range(num_fem):
	o.write('female%s female%s female%s 2\n' % (ix, ix, ix))
o.close()

legendfile = '%sfemalehaps.legend' % outdir
o = open(legendfile, 'w')
o.write('id position a0 a1\n')
for ix, pos in enumerate(sorted(snps.keys())):
	o.write('SNP%s %s %s %s\n' % (ix, pos, snps[pos]['ref'], snps[pos]['alt']))
o.close()

hapfile = '%sfemalehaps.hap' % outdir
o = open(hapfile, 'w')
for pos in sorted(snps.keys()):
	o.write(' '.join(snps[pos]['gen']) + '\n')
o.close()
	
