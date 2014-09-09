import re
import random

file = '/mnt/gluster/home/sonal.singhal1/Darwin/g_fortis/Geospiza_fortis.fasta'

hets = { 'M': ['A', 'C'], 'R': ['A', 'G'], 'W': ['A', 'T'], 'S': ['C', 'G'], 'Y': ['C', 'T'], 'K': ['G', 'T'] }

dir = '/mnt/gluster/home/sonal.singhal1/gene_trees/chromosomes/'

f_in = open(file, 'r')
chromosome = ''
for l in f_in:
	if re.match('>', l):
		if len(chromosome) > 0:
			haplo0 = list(chromosome)
			haplo1 = list(chromosome)
			for ix, bp in enumerate(haplo0):
				if bp in hets:
					if random.random() > 0.5:
						haplo0[ix] = hets[bp][0]
						haplo1[ix] = hets[bp][1]
					else:
						haplo0[ix] = hets[bp][1]
                                                haplo1[ix] = hets[bp][0]
			f_out.write('>haplo0\n')
			for i in xrange(0, len(haplo0), 60):
	                	f_out.write(''.join(haplo0[i:i+60]) + '\n')
			f_out.write('>haplo1\n')
                        for i in xrange(0, len(haplo1), 60):
                                f_out.write(''.join(haplo1[i:i+60]) + '\n')
			haplo0 = ''
			haplo1 = ''
			chromosome = ''
			f_out.close()
		chr = re.match('>(\S+)', l).group(1)
		file = '%sGeospiza_%s_haplotypes.fasta' % (dir, chr)
		f_out = open(file, 'w')
	else:
		chromosome += l.rstrip().upper()
f_in.close()
