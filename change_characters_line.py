import re

ref_genome = '/mnt/gluster/home/sonal.singhal1/Darwin/g_fortis/Geospiza_singleline.fa'
out = '/mnt/gluster/home/sonal.singhal1/Darwin/g_fortis/geo60.fa'

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
		left += l
		if len(left) >= 60:
			outfile.write("%s\n" % left[0:60])
			left = left[60:]
outfile.write("%s\n" % left)
infile.close()
outfile.close()
			
