import re

chrs = ['chr1', 'chr1A', 'chr3', 'chr4', 'chr4A', 'chr5', 'chr6',
        'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 
        'chr15', 'chr2']

compare = {'ZF-Geo': ['zebrafinch', 'geospiza'], 'ZF-Fic': ['zebrafinch', 'ficedula']}

legal = ['A', 'T', 'C', 'G']
out = '/mnt/gluster/home/sonal.singhal1/MSA/divergences.csv'
o = open(out, 'w')
o.write('comparison,chr,num_diff,num_sites\n')

for chr in chrs:
	for comparison in compare:
		sp1 = compare[comparison][0]
		sp2 = compare[comparison][1]

		maf = '/mnt/gluster/home/sonal.singhal1/MSA/%s/%s.%s.sing.maf' % (chr, sp1, sp2)

		num_sites = 0
		num_diff = 0

		f = open(maf, 'r')
		for l in f:
			if re.match('a score=\d+', l):
				seq1 = re.split('\s+', f.next().rstrip())[6]
				seq2 = re.split('\s+', f.next().rstrip())[6]

				for bp1, bp2 in zip(seq1, seq2):
					# okay not to do soft-masking support
					# tba doesn't handle soft-masked sequences well
					if bp1 in legal and bp2 in legal:
						num_sites += 1
						if bp1 != bp2:
							num_diff += 1
		o.write('%s,%s,%s,%s\n' % (comparison, chr, num_diff, num_sites))
