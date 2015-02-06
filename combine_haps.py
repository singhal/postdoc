import re
import glob

chr = 'chr2'
files = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/finch21/%s*haps' % chr)
out = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/finch21/%s_haplotypes.haps'  % chr

sites = {}

for file in files:
	f = open(file, 'r')
	for l in f:
		d = re.split('\s+', l)
		if int(d[2]) not in sites:
			sites[ int(d[2]) ] = l
	f.close()

o = open(out, 'w')
ordered_sites = sorted(sites.keys())
for site in ordered_sites:
	o.write(sites[site])
o.close()
