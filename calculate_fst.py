import re
import gzip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--chr", help="chromosome for which to run analysis")
args = parser.parse_args()
chr = args.chr

vcf = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/for_shapeit/gatk.ug.ltf.%s.allfilters.recoded_biallelicSNPs.vcf.gz'  % chr
if chr == 'chrZ':
	vcf = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/for_shapeit/gatk.ug.ltf.chrZ.allfilters.recodedsex.recoded_biallelicSNPs.males.vcf.gz'
out = '/mnt/gluster/home/sonal.singhal1/LTF/analysis/fst/%s.Wrights_Fst.csv.gz' % chr

o = gzip.open(out, 'w')
o.write('chr,pos,Fst\n')
def get_p(pop):
	num_p = 0
	for gen in pop:
		if gen == '0/0':
			num_p += 2
		elif gen in ['1/0', '0/1']:
			num_p += 1
	af_p = num_p / float(len(pop) * 2)
	return af_p

f = gzip.open(vcf, 'r')
for l in f:
	if not re.search('#', l):
		d = re.split('\t', l.rstrip())
		genos = []
		for i in d[9:]:
			geno = re.search('^(\S/\S)', i).group(1)
			genos.append(geno)
		
		# account for chromosome Z
		if chr != 'chrZ':
			pop1 = genos[0:10]
			pop2 = genos[10:20]
		else:
			pop1 = genos[0:5]
			pop2 = genos[5:12]

		# account for missing data
		pop1 = [x for x in pop1 if not re.search('\.', x)]
		pop2 = [x for x in pop2 if not re.search('\.', x)]
		genos = pop1 + pop2

		if len(pop1) > 0 and len(pop2) > 0:
			pop1_p = get_p(pop1)
			pop2_p = get_p(pop2)
			tot_p = get_p(genos)

			hs_1 = 2 * pop1_p * (1 - pop1_p)
			hs_2 = 2 * pop2_p * (1 - pop2_p)
			ht = 2 * tot_p * (1 - tot_p)

			if ht != 0:
				fst = (ht - ((hs_1 + hs_2) / 2.0))/ ht
				o.write('%s,%s,%s\n' % (chr, d[1], fst))
			else:
				o.write('%s,%s,%s\n' % (chr, d[1], 'nan'))
		else:
			o.write('%s,%s,%s\n' % (chr, d[1], 'nan'))
o.close()
