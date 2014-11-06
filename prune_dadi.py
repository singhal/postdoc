import re

file = 'all_species.dadi.txt'
outfile = 'all_species.dadi_pruned.txt'
out = open(outfile, 'w')

f = open(file, 'r')
out.write(f.next().rstrip() + '\n')
pos = 0
chr = 'chrLGE22'
for l in f:
	l = l.rstrip()
	d = re.split('\t', l)
	if int(d[11]) - pos > 3:
		out.write(l + '\n')
		pos = int(d[11])
	if d[10] != chr:
		pos = 0
		chr = d[10]
out.close()
