import re

ref_genome = '/mnt/lustre/home/sonal.singhal1/reference/taeGut1.bamorder.fasta'
out = '/mnt/lustre/home/sonal.singhal1/reference/taeGut1_60.bamorder.fasta'

infile = open(ref_genome, 'r')
outfile = open(out, 'w')

left = ''
for l in infile:
	l = l.strip()
	if re.match(">", l):
		outfile.write("%s\n" % left)
		left = ''
		outfile.write("%s\n" % l)
	else:
		l = l.upper()
		left += l
		if len(left) >= 60:
			outfile.write("%s\n" % left[0:60])
			left = left[60:]
outfile.write("%s\n" % left)
infile.close()
outfile.close()
			
