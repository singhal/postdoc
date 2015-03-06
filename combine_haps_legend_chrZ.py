import re
import gzip

sp = 'LTF'

if sp == 'ZF':
	legend = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/all_vcf/for_shapeit/chrZ.hap.gz'
	haps = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/finch19/chrZ_haplotypes.haps'
if sp == 'LTF':
	legend = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/for_shapeit/chrZ.hap.gz'
	haps = '/mnt/gluster/home/sonal.singhal1/LTF/phasing/PIR_approach/chrZ_haplotypes.haps' 

haps_data = []
f = gzip.open(legend, 'r')
for l in f:
	d = re.split('\s', l.rstrip())
	haps_data.append(d[0::2])
f.close()

out = haps + '2'
o = open(out, 'w')
f = open(haps, 'r')
for hap, l in zip(haps_data, f):
	l = l.rstrip()
	o.write(l + ' '  + ' '.join(hap) + '\n')
f.close()
o.close()

	
