import glob
import re

files = glob.glob('/mnt/gluster/home/sonal.singhal1/gene_trees2/individual_trees/fasta_files/*aln')

comparisons = {'ZF-LTF': ['ZF_26792b', 'LTFh_73942b'], 'ZF-DBF': ['ZF_26792b', 'DBFa'], 'ZF-Ficedula': ['ZF_26792b', 'ficedula'], 'ZF-Geospiza': ['ZF_26792b', 'gfortis']}
out = '/mnt/gluster/home/sonal.singhal1/gene_trees2/divergences.csv'
o = open(out, 'w')
o.write('gene,comparison,num_sites,num_diff\n')

for file in files:
	name = re.search('(chr\S+)\.fasta', file).group(1)
	
	seq = {}
	id = ''
	
	f = open(file, 'r')
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l).group(1)
		else:
			if id not in seq:
				seq[id] = ''
			seq[id] += l.rstrip()
	f.close()

	for compare in comparisons:
		sp1 = comparisons[compare][0]
		sp2 = comparisons[compare][1]

		denom = 0
		num_diff = 0

		for a, b in zip(seq[sp1], seq[sp2]):
			if a in ['A', 'T', 'C', 'G'] and b in ['A', 'T', 'C', 'G']:
				denom += 1
				if a != b:
					num_diff += 1

		o.write('%s,%s,%s,%s\n' % (name, compare, num_diff, denom))
o.close()
