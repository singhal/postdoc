import re
import glob

files = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/mendelian_errors/plink_results/*.mendel')

mendel = {}

for file in files:
	chr = re.search('(chr[A-Z|0-9]+)', file).group(1)
	mendel[chr] = {'snp': 0, 'indel': 0}
	
	f = open(file, 'r')
	head = f.next()
	for l in f:
		d = re.split('\s\s+', l.rstrip())
		mut = d[-1]
		type = 'indel'
		if re.search('[A|T|G|C]', mut):
			type = 'snp'
		mendel[chr][type] +=1
	f.close()

print 'chr,snp,indel'
for chr in mendel:
	print '%s,%s,%s' % (chr, mendel[chr]['snp'], mendel[chr]['indel']) 
