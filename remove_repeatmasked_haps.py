import re

chrs = ['chr3', 'chr13', 'chr19', 'chr20', 'chr22', 'chr23', 'chr28']
repeat = '/mnt/gluster/home/sonal.singhal1/reference/taeGut1.repeatMaskerBlast.repeatLibrary20140131.out'

repeats = {}
f = open(repeat, 'r')
for l in f:
        if re.search('^\s+\d+', l):
                d = re.split('\s+', l.rstrip())
                if d[5] in chrs:
                        if d[5] not in repeats:
                                repeats[d[5]] = {}
                        start = int(d[6])
                        end = int(d[7]) + 1

                        for pos in range(start, end):
                                repeats[d[5]][pos] = 1
f.close()

for chr in chrs:
	old = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/old_results/%s_haplotypes.haps' % chr
	new = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/results/%s_haplotypes.haps' % chr

	f = open(old, 'r')
	o = open(new, 'w')

	for l in f:
		d = re.split('\s+', l)
		if int(d[2]) not in repeats[chr]:
			o.write(l)
	o.close()
	f.close()
