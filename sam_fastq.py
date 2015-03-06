import re

f = open('Geospiza_fortis.unmapped.sorted.sam', 'r')
o1 = open('one.fastq', 'w')
o2 = open('two.fastq', 'w')
u = open('unpaired.fastq', 'w')

def print_line(d, out):
	out.write('@' + d[0] + '\n')
	out.write(d[9] + '\n')
	out.write('+\n')
	out.write(d[10] + '\n')

id1 = ''
line1 = []
for l in f:
	d = re.split('\s+', l.rstrip())
	id = re.sub('_\d', '',  d[0])
	if id1:
		if id == id1:
			print_line(line1, o1)
			print_line(d, o2)
			id1 = ''
			line1 = []
		else:
			print_line(line1, u)
			id1 = id
			line1 = d
	else:
		id1 = id
		line1 = d
o1.close()
f.close()
o2.close()
u.close()	
