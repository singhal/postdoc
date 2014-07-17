from itertools import izip
import gzip
import re
import sys

def get_vcf(vcffile):
	file = gzip.open(vcffile, 'r')
	var = dict()
	for l in file:
		if not re.match("#", l):
			d = re.split('\t', l)

			if d[0] not in var:
				var[d[0]] = dict()
			alleles = [d[3]] + re.split(",", d[4])
			allele_counts = dict()
			for i in range(len(alleles)):			
				allele_counts[alleles[i]] = len(re.findall('%s\/' % i, l)) + len(re.findall('\/%s' % i, l))
			tot_n = sum(allele_counts.values())
			for i in allele_counts:
				allele_counts[i] = '%.3f' % (allele_counts[i] / float(tot_n))
			var[d[0]][int(d[1])] = allele_counts
	file.close()
	return var	


vcf3 = '/mnt/lustre/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.chr1.filtered.coverage.vqsr.vcf.test.gz'
vcf2 = '/mnt/lustre/home/sonal.singhal1/ZF/after_vqsr/by_chr/gatk.ug.zf.chr1.filtered.coverage.vqsr.vcf.test.gz'
vcf1 = '/mnt/lustre/home/sonal.singhal1/DBF/after_vqsr/by_chr/gatk.ug.dbf.chr1.filtered.coverage.vqsr.vcf.test.gz'

# importantly, all these genomes are interleaved at 60chr long
# this prevents me from having to read the entire genome into memory
# should be much faster than the access methods too
#
# masked genome files for each species; genome1 and genome 2 are outgroup
genome1 = '/mnt/lustre/home/sonal.singhal1/DBF/masked_genome/DBF.masked_genome.fa'
genome2 = '/mnt/lustre/home/sonal.singhal1/ZF/masked_genome/ZF.masked_genome.fa'
# genome 3 is the species for which we are generating the mutation matrix
genome3 = '/mnt/lustre/home/sonal.singhal1/LTF/masked_genome/LTF.masked_genome.fa'

# the genome that tells the reference site at each position
# need this to get site count of the entire genome
ref_genome = '/mnt/lustre/home/sonal.singhal1/reference/taeGut1_60.bamorder.fasta'

var1 = get_vcf(vcf1)
var2 = get_vcf(vcf2)
var3 = get_vcf(vcf3)

f1 = open(genome1, 'r')
f2 = open(genome2, 'r')
f3 = open(genome3, 'r')
ref = open(ref_genome, 'r')

line_count = 0
for l1, l2, l3, refline in izip(f1, f2, f3, ref):
	l1 = l1.strip()
	l2 = l2.strip()
	l3 = l3.strip()
	refline = refline.strip()

	if re.match(">(\S+)", l1):
		chr = re.match(">(\S+)", l1).group(1)
		# vcf = dir + 'gatk.ug.ltf.%s.filtered.coverage.biallelicSNPs.vqsr.vcf.gz' % chr
		
		# need to empty out the variation dictionary
		# var = dict()
		#  need to reset the line count
		line_count = 0
	
		# vcfopen = gzip.open(vcf, 'r')
		# for l in vcfopen:
		#	if not re.match("#", l):
		#		d = re.split('\t', l)
		#		var[int(d[1])] = {'ref': d[3], 'alt': d[4]}
		# vcfopen.close()
	else:
		for ix, (a, b, c, d) in enumerate(izip(l1,l2,l3,refline)):	
			a = int(a)
			b = int(b)
			c = int(c)

			# get the current genome position	
			# need to add 1 because python indexes at 0; genomes are at index 1
			pos = line_count * 60 + ix + 1
			if c == 1:
				print var3[chr][pos]
				if a == 1:
					print var1[chr][pos]
				else:
					if a < 4:
						print "MONO"
					else:
						print "NONE"
				if b == 1:
					print var2[chr][pos]
				else:
					if b < 4:
						print "MONO"
					else:
						print "NONE"
				print "******"
				

			if pos > 500000:
				sys.exit()

		line_count += 1

f1.close()
f2.close()
f3.close()
ref.close()
out.close()
