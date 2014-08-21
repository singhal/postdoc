import re
import glob
import gzip

dir = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/'
files = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/*vcf.gz')

for file in files:
	chr = re.search('(chr[a-z|A-Z|0-9]+)', file).group(1)
	out = '%sgatk.ug.all_zf.%s.coverage.filtered.vqsr.vcf.gz' % (dir, chr)
	
	out_f = gzip.open(out, 'w')
	in_f = gzip.open(file, 'r')
	
	for l in in_f:
		if re.search('^#', l):
			out_f.write(l)
		else:
			d = re.split('\s+', l)
			if d[4] !=  '.':
				if re.search('PASS', d[6]):
					out_f.write(l)
	out_f.close()
	in_f.close()
