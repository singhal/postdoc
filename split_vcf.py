import re

file = 'Geospiza_fortis.vcf'
f = open(file, 'r')

chr = ''

for l in f:
	if re.search('^chr', l):
		d = re.split('\t', l)
		if d[0] != chr:
			chr = d[0]
			try:
				o.close()
			except:
				pass
			out = file.replace('.vcf', '_%s.vcf' % chr)
			o = open(out, 'w')
		o.write(l)
f.close()
