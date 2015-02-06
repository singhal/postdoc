import re
import pandas as pd

# get gene names from groundfinch
gf = {}
gf_gff = '/mnt/gluster/home/sonal.singhal1/PAR/groundfinch/Geospiza_fortis.gene.gff'
f = open(gf_gff, 'r')
for l in f:
	d = re.split('\t', l)
	if d[0] == 'scaffold143':
		if re.search('Function', d[8]):
			genename = re.search('Function=\"(.*)\"', d[8]).group(1)
			geneid = re.search('Source=([A-Z|0-9]+)', d[8]).group(1)
			gf[geneid] = genename
f.close()

# get gene names from flycatcher
fc = {}
fc_names = '/mnt/gluster/home/sonal.singhal1/PAR/flycatcher/flycatcher_PAR_genes.csv'
d = pd.read_csv(fc_names)
for name in d.GeneSymbol:
	fc[name] = 1

# get gene names from zebrafinch
zf_gff = '/mnt/gluster/home/sonal.singhal1/reference/Taeniopygia_guttata.gff'
f = open(zf_gff, 'r')
zf = {}
for l in f:
	d = re.split('\s+', l)
	if d[0] == 'Z_random':
		if int(d[4]) < 500000:
			if re.search('ID', l):
				geneid = re.search('ID=(.*);', l).group(1)
				zf[geneid] = 1
f.close()

zf_names = '/mnt/gluster/home/sonal.singhal1/reference/Ensembl_data_ZF_genes.csv'
f = open(zf_names, 'r')
for l in f:
	d = re.split(',', l.rstrip())
	if d[3] in zf:
		names = [d[1], d[2], d[4]]
		names = filter(lambda x: not re.match('LOC\d+', x), names)
		names = filter(lambda x: not re.match('^$', x), names)
		if len(names) < 1:
			names = ['NA']
		zf[d[3]] = names[0]
f.close()

for contig in fc:
	print 'fc\t%s' % contig
for id, name in gf.items():
	print 'gf\t%s\t%s' % (name, id)
for id, name in zf.items():
	print 'zf\t%s\t%s' % (name, id)

