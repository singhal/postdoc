import re
import gzip

dir = '/mnt/gluster/home/sonal.singhal1/inversions/'
file = '%sEstrildidFinch_INVY.vcf.gz' % dir
# file = '%stest.vcf.gz' % dir

def af(genos):
	genos = [x for x in genos if x != '.']

	af = {}
	for allele in ['0', '1']:
		if len(genos) > 10:
			af[allele] = genos.count(allele) / float(len(genos))
		
	return af, len(genos)

f = gzip.open(file, 'r')
for l in f:
	if not re.match('^#', l):
		d = re.split('\t', l.rstrip())
		if d[6] == 'PASS':
			genos = []
			for geno in d[9:]:
				geno = re.search('^(\S\S\S)', geno).group(1)
				geno = re.split('/', geno)
				genos += geno
			zf = genos[0:38]
			ltf = genos[38:]

			af_zf, zfn = af(zf)
			af_ltf, ltfn = af(ltf)

			fixed = False
			for allele in ['0', '1']:
				if allele in af_zf and allele in af_ltf:
					diff = abs(af_zf[allele] - af_ltf[allele])
					if diff > 0.95:
						fixed = True

			pe = int(re.search('PE=(\d+)', d[7]).group(1))
			mapq = int(re.search('MAPQ=(\d+)', d[7]).group(1))
			
			
			if fixed and pe > 50 and mapq > 30:
				print '%s %.2f %s %s' % (d[0], diff, zfn, ltfn)
