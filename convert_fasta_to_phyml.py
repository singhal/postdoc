import glob
import re

files = glob.glob('/mnt/gluster/home/sonal.singhal1/gene_trees/individual_trees/phyml_files/*phyml')

for file in files:
	'''
	f_in = open(file, 'r')
	out = file.replace('fasta', 'phyml')
	seq = {}
	id = ''
	for l in f_in:
		if re.match('>', l):
			id = re.match('>(\S+)', l).group(1)
		else:
			if id not in seq:
				seq[id] = ''
			seq[id] += l.rstrip()
	f_in.close()

	out_f = open(out, 'w')
	out_f.write('%s %s\n' % (len(seq), len(seq[seq.keys()[0]])))
	for id, sequence in seq.items():
		out_f.write('%s %s\n' % (id, sequence))
	out_f.close()
	'''
	print '%s %s' % ('~/bin/mraic.pl', file)

