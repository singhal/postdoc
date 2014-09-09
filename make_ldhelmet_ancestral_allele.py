import re
import glob

chrs = [ 'chr1', 'chr1A', 'chr1B', 'chr2', 'chr3',  'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
         'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
         'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28',
         'chrLG2', 'chrLG5', 'chrLGE22', 'chrZ' ]

nuc = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

for chr in chrs:
	aa_file = '/mnt/gluster/home/sonal.singhal1/ZF/ancestral_allele/ancestral_allele.%s.csv' % chr
	sites_file = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/%s_sites.csv' % chr
	out_file = '/mnt/gluster/home/sonal.singhal1/ZF/ancestral_allele/ancestral_allele.%s.ldhelmet.txt' % chr

	anc_info = {}
	
	aa_f = open(aa_file, 'r')
	header = aa_f.next()
	for l in aa_f:
		d = re.split(',', l.rstrip())
		site = int(d[1])

		if d[-1] == 'N':
			anc_info[site] = [0.30, 0.20, 0.20, 0.30]
		else:
			anc_info[site] = [0.02, 0.02, 0.02, 0.02]
			anc_info[site][nuc[d[-1]]] = 0.94
	aa_f.close()

	out_f = open(out_file, 'w')
	sites_f = open(sites_file, 'r')
	for l in sites_f:
		d = re.split(',', l.rstrip())
		d[1] = int(d[1])
		if d[1] in anc_info:
			out_f.write('%s %s\n' % (d[1] - 1, ' '.join('%.2f' % i for i in anc_info[int(d[1])])))
		else:
			out_f.write('%s %s\n' % (d[1] - 1, '0.29 0.21 0.21 0.29'))
			print 'Undefined ancestral sequence: %s %s' % (chr, d[1])
	out_f.close()
	sites_f.close()
