import re
#from collections import Counter

file = '/mnt/gluster/home/sonal.singhal1/reference/taeGut1.bamorder.fasta'

f_open = open(file, 'r')
sites = {'A': 0, 'G': 0, 'N': 0, 'T': 0, 'C': 0}

for l in f_open:
	if not re.search('^>', l):
		l = list(l.rstrip())
		for base in sites:
			sites[base] += l.count(base)

for base, count in sites.items():
	print '%s %s' % (base, count)


