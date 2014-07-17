from itertools import izip
import re
import gzip

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

# the drive that has the vcf files with the variant info
dir = '/mnt/lustre/home/sonal.singhal1/LTF/after_vqsr/by_chr/'
# file that will have all the info to generate the mutation matrix
outfile = dir + 'mutation_matrix_info.txt'
out = open(outfile, 'w')

f1 = open(genome1, 'r')
f2 = open(genome2, 'r')
f3 = open(genome3, 'r')
ref = open(ref_genome, 'r')

# initialize the sites that store all the info
site_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
mut_matrix = {}
for bp1 in ['A', 'T', 'C', 'G']:
	mut_matrix[bp1] = dict()
	for bp2 in ['A', 'T', 'C', 'G']:
		mut_matrix[bp1][bp2] = 0

line_count = 0
for l1, l2, l3, refline in izip(f1, f2, f3, ref):
	l1 = l1.strip()
	l2 = l2.strip()
	l3 = l3.strip()
	refline = refline.strip()

	if re.match(">(\S+)", l1):
		chr = re.match(">(\S+)", l1).group(1)
		vcf = dir + 'gatk.ug.ltf.%s.filtered.coverage.biallelicSNPs.vqsr.vcf.gz' % chr
		
		# need to empty out the variation dictionary
		var = dict()
		# need to reset the line count
		line_count = 0
	
		vcfopen = gzip.open(vcf, 'r')
		for l in vcfopen:
			if not re.match("#", l):
				d = re.split('\t', l)
				var[int(d[1])] = {'ref': d[3], 'alt': d[4]}
		vcfopen.close()
	else:
		for ix, (a, b, c, d) in enumerate(izip(l1,l2,l3,refline)):	
			a = int(a)
			b = int(b)
			c = int(c)

			# get the current genome position	
			# need to add 1 because python indexes at 0; genomes are at index 1
			pos = line_count * 60 + ix + 1
		
			if a < 4 and b < 4 and c < 4:
				# only count sites that are at sufficient coverage in all three genomes
				site_counts[d] += 1
				# only count sites that are variable in the ref genome
				if a != 1 and b != 1:
					if c == 1:
						# not all variable positions are in the var file
						# because var file only includes biallelic; that's all LDhelmet can take
						if pos in var:
							mut_matrix[ var[pos]['ref'] ][ var[pos]['alt'] ] += 1

		line_count += 1

for bp in site_counts:
	out.write("sitecounts\t%s\t%s\n" % (bp, site_counts[bp]))
for bp1 in mut_matrix:
	for bp2 in mut_matrix[bp1]:
		out.write("mutmatrix\t%s_%s\t%s\n" % (bp1, bp2, mut_matrix[bp1][bp2]))
f1.close()
f2.close()
f3.close()
ref.close()
out.close()
