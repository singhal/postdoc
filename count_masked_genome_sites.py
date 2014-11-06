import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--genome", help="genome for which you want to evaluate sites")

args = parser.parse_args()
genome = args.genome

out = genome + '.site_counts.csv'
total = {}
chr = ''
counts = {}

o = open(out, 'w')
o.write('chr,site_type,count\n')
gen = open(genome, 'r')
for l in gen:
	l = l.rstrip()
	if re.search('>', l):
		if chr:
			for ix, count in counts.items():
				o.write('%s,%s,%s\n' % (chr, ix, count))
				if ix not in total:
					total[ix] = 0
				total[ix] += count
			counts = {}
		chr = re.search('>(\S+)', l).group(1)
	else:
		for bp in l:
			if bp not in counts:
				counts[bp] = 0
			counts[bp] += 1
gen.close()
for ix, count in total.items():
	o.write('%s,%s,%s\n' % ('total', ix, count))
o.close()			
