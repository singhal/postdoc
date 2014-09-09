import re

file = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/PSMC/26462.fasta'

f_open = open(file, 'r')
sites = {'A': 0, 'G': 0, 'T': 0, 'C': 0, 'M': 0, 'R': 0, 'W': 0, 'S': 0, 'Y': 0, 'K':0}

for l in f_open:
	if not re.search('^>', l):
		l = list(l.replace('N', '').rstrip())
		for base in sites:
			sites[base] += l.count(base)

for base, count in sites.items():
	print '%s %s' % (base, count)


