import re
import argparse
import gzip

def get_var(vcf, chr, length):
	seq = ['N'] * length

	ambig = {	'A': {'C': 'M', 'G': 'R', 'T': 'W', 'N': 'N'},
			'G': {'A': 'R', 'C': 'S', 'T': 'K', 'N': 'N'},
			'C': {'A': 'M', 'G': 'S', 'T': 'Y', 'N': 'N'},
			'T': {'A': 'W', 'C': 'Y', 'G': 'K', 'N': 'N'},
			'N': {'A': 'N', 'C': 'N', 'G': 'N', 'T': 'N', 'N': 'N'}
		}

	f = gzip.open(vcf, 'r')
	for l in f:
		if not re.search('^#', l):
			d = re.split('\t', l.rstrip())
			depth = int(re.search('DP=(\d+)', l).group(1))
			if depth >= 3:
				pos = int(d[1]) - 1
				if d[4] == '.':
					seq[pos] = d[3]
				else:
					alleles = [d[3]] + re.split(',', d[4])
					genos = {}
					num = -1
					for ix, allele in enumerate(alleles):
						for sec_allele in alleles[0:ix]:
							num = num + 1
							genos[num] = ambig[allele][sec_allele]
						num = num + 1	
						genos[num] = allele			

					lk = re.split(',', d[-1])

					geno = genos[lk.index(min(lk))]
					seq[pos] = geno

	f.close()
	return seq						

				
def print_seq(chr, seq, out):			
	o = open(out, 'w')
	o.write('>%s\n' % chr)
        for i in xrange(0, len(seq), 60):
                o.write(''.join(seq[i:i+60]) + '\n')
        o.close()


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--chr", help="chromosome for which to run analysis")
	parser.add_argument("--sp", help="species for which to run analysis")
	args = parser.parse_args()
	chr = args.chr
	sp = args.sp

	if sp == 'gmag':
		dir = 'g_magnirostris'
		name = 'Geospiza_magnirostris'
	elif sp == 'gfor':
		dir = 'g_fortis'
		name = 'Geospiza_fortis'
	
	chr_lengths = { 'chr10': 20806668, 'chr11': 21403021, 'chr12': 21576510, 'chr13': 16962381,
                'chr14': 16419078, 'chr15': 14428146, 'chr16': 9909, 'chr17': 11648728,
                'chr18': 11201131, 'chr19': 11587733, 'chr1A': 73657157, 'chr1B': 1083483,
                'chr1': 118548696, 'chr20': 15652063, 'chr21': 5979137, 'chr22': 3370227,
                'chr23': 6196912, 'chr24': 8021379, 'chr25': 1275379, 'chr26': 4907541,
                'chr27': 4618897, 'chr28': 4963201, 'chr2': 156412533, 'chr3': 112617285,
                'chr4A': 20704505, 'chr4': 69780378, 'chr5': 62374962, 'chr6': 36305782,
                'chr7': 39844632, 'chr8': 27993427, 'chr9': 27241186, 'chrLG2': 109741,
                'chrLG5': 16416, 'chrLGE22': 883365, 'chrZ': 72861351, 'chrZ_random': 2969867}

	length = chr_lengths[chr]
	vcf = '/mnt/gluster/home/sonal.singhal1/ficedula/Ficedula_albicollis.%s.vcf.gz' % chr
	out = '/mnt/gluster/home/sonal.singhal1/ficedula/Ficedula_albicollis.%s.fa' % chr
	#vcf = '/mnt/gluster/home/sonal.singhal1/Darwin/%s/vcf/%s_%s.vcf.gz' % (dir, name, chr)
	#out = '/mnt/gluster/home/sonal.singhal1/Darwin/%s/%s_%s.fa' % (dir, name, chr)

	seq = get_var(vcf, chr, length)
	print_seq(chr, seq, out)

if __name__ == "__main__":
    main()
