import re
import gzip
import glob
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--sp", help="species for which to run analysis")
args = parser.parse_args()
sp = args.sp

chr_lengths = {	'chr10': 20806668, 'chr11': 21403021, 'chr12': 21576510, 'chr13': 16962381,
                'chr14': 16419078, 'chr15': 14428146, 'chr16': 9909, 'chr17': 11648728,
                'chr18': 11201131, 'chr19': 11587733, 'chr1A': 73657157, 'chr1B': 1083483,
                'chr1': 118548696, 'chr20': 15652063, 'chr21': 5979137, 'chr22': 3370227,
                'chr23': 6196912, 'chr24': 8021379, 'chr25': 1275379, 'chr26': 4907541,
                'chr27': 4618897, 'chr28': 4963201, 'chr2': 156412533, 'chr3': 112617285,
                'chr4A': 20704505, 'chr4': 69780378, 'chr5': 62374962, 'chr6': 36305782,
                'chr7': 39844632, 'chr8': 27993427, 'chr9': 27241186, 'chrLG2': 109741,
                'chrLG5': 16416, 'chrLGE22': 883365, 'chrZ': 72861351}
genome = '/mnt/gluster/home/sonal.singhal1/%s/masked_genome/%s.masked_genome.repeat_masked.fa' % (sp, sp)
window = 50000

if sp == 'ZF':
	nind = {}
	for chr in chr_lengths:
		nind[chr] = 38
	nind['chrZ'] = 28
	vcfs = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/*phased*')
if sp == 'LTF':
	nind = {}
        for chr in chr_lengths:
                nind[chr] = 40
        nind['chrZ'] = 32
	vcfs = glob.glob('/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/*phased*')
theta = {}
outfile = '/mnt/gluster/home/sonal.singhal1/%s/analysis/pop_gen/wattersons_theta.csv' % sp

for chr, length in chr_lengths.items():
	theta[chr] = {}
	for ix, start in enumerate(range(1, length, window)):
		end = start + window
		if end > length:
			end = length
		
		seq = subprocess.Popen('samtools faidx %s %s:%s-%s' % (genome, chr, start, end), shell=True, stdout=subprocess.PIPE)
		n_count = 0
		for l in seq.stdout:
			if not re.match('>', l):
				n_count += len(filter(lambda x: int(x) < 4, list(l.rstrip())))

		theta[chr][ix] = {'start': start, 'end': end, 'ss': 0, 'seq_length': n_count}


for gzvcf in vcfs:
	f = gzip.open(gzvcf, 'r')
	chr = re.search('(chr[A-Z|0-9]+)', gzvcf).group(1)
	for l in f:
		if not re.search('^#', l):
			d = re.split('\t', l)
			# look for non-indels
			alleles = [d[3]] + re.split(',', d[4])
			indel = False
			for allele in alleles:
				if len(allele) > 1:
					indel = True
			if not indel:
				genos = []
				for geno in d[9:]:
					geno = re.search('^([^:]+)', geno).group(1)
					genos += re.split('[|/]', geno)
				num_alleles = 0
				for ix, allele in enumerate(alleles):
					if genos.count(str(ix)) > 0:
						num_alleles += 1
				if num_alleles > 1:
					# look for variable sites only
					# identifies the chunk into which the variant goes
					ix = int(int(d[1])/float(window))
					theta[chr][ix]['ss'] += 1
	f.close()

out = open(outfile, 'w')
out.write('chr,index,start,end,seq_length,watterson_theta\n')
for chr in theta:

	harmonic = 0
	# to account for different nind in different chromosomes
	for i in range(1, nind[chr]):
		harmonic += 1/float(i)

	for ix in sorted(theta[chr].keys()):
		if theta[chr][ix]['seq_length'] > 0:
			thetaw = theta[chr][ix]['ss'] / float(theta[chr][ix]['seq_length']) / harmonic
			out.write('%s,%s,%s,%s,%s,%.4f\n' % (chr, ix, theta[chr][ix]['start'], theta[chr][ix]['end'], theta[chr][ix]['seq_length'], thetaw))
		else:
			out.write('%s,%s,%s,%s,0,NA\n' % (chr, ix, theta[chr][ix]['start'], theta[chr][ix]['end']))
out.close()
