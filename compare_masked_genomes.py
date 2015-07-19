import re
from itertools import izip

mgen1 = '/mnt/gluster/home/sonal.singhal1/DBF/masked_genome/DBF.masked_genome.fa'
mgen2 = '/mnt/gluster/home/sonal.singhal1/LTF/masked_genome/LTF.masked_genome.repeat_masked.fa'
mgen3 = '/mnt/gluster/home/sonal.singhal1/ZF/masked_genome/ZF.masked_genome.repeat_masked.switch_masked.fa'

out = '/mnt/gluster/home/sonal.singhal1/reference/all.masked_genome.fa'
f_out = open(out, 'w')

m_f1 = open(mgen1, 'r')
m_f2 = open(mgen2, 'r')
m_f3 = open(mgen3, 'r')

for l1, l2, l3 in izip(m_f1, m_f2, m_f3):
	if re.match('>', l1):
		f_out.write(l1)
	else:
		l = []
		l1 = [int(x) for x in list(l1.rstrip())]
		l2 = [int(x) for x in list(l2.rstrip())]
		l3 = [int(x) for x in list(l3.rstrip())]
	
		for bp1, bp2, bp3 in izip(l1, l2, l3):
			num_right = 0
			if bp1 < 4:
				num_right += 1
			if bp2 < 4:
				num_right += 1
			if bp3 < 4:
				num_right += 1
			if num_right >= 3:
				l.append('0')
			else:
				l.append('8')
		f_out.write(''.join(l) + '\n')
m_f1.close()
m_f2.close()
m_f3.close()
f_out.close()
