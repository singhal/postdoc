import gzip
import re
import glob

files = glob.glob('/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/for_shapeit/*vcf*gz')

out = open('LTF.shapeitVCFs_counts.csv', 'w')
out.write('file,count\n')

for file in files:
	f = gzip.open(file, 'r')
	ix = 0
	for l in f:
		if not re.match('^#', l):
			ix += 1
	f.close()
	out.write('%s,%s\n' % (file, ix))
out.close()
