import glob
import gzip
import re

dir = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/all_vcf/'
files = glob.glob(dir + '*filt*vqsr2*')

out = '%sgatk.ug.zf.all_chrs.masked.filtered.vqsr2.vcf.gz' % dir
o = open(out, 'w')

for ix, file in enumerate(files):
	print file
	if ix == 0:
		header = True
	else:
		header = False

	f = gzip.open(file, 'r')
	for l in f:
		if re.search('^#', l):
			if header:
				o.write(l)
		else:
			o.write(l)
	f.close()
o.close()
print 'done!'
