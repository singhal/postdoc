import re
import gzip
import glob
import subprocess

chr_lengths = {	'chr10': 20806668, 'chr11': 21403021, 'chr12': 21576510, 'chr13': 16962381,
                'chr14': 16419078, 'chr15': 14428146, 'chr16': 9909, 'chr17': 11648728,
                'chr18': 11201131, 'chr19': 11587733, 'chr1A': 73657157, 'chr1B': 1083483,
                'chr1': 118548696, 'chr20': 15652063, 'chr21': 5979137, 'chr22': 3370227,
                'chr23': 6196912, 'chr24': 8021379, 'chr25': 1275379, 'chr26': 4907541,
                'chr27': 4618897, 'chr28': 4963201, 'chr2': 156412533, 'chr3': 112617285,
                'chr4A': 20704505, 'chr4': 69780378, 'chr5': 62374962, 'chr6': 36305782,
                'chr7': 39844632, 'chr8': 27993427, 'chr9': 27241186, 'chrLG2': 109741,
                'chrLG5': 16416, 'chrLGE22': 883365, 'chrZ': 72861351}
genome = '/mnt/gluster/home/sonal.singhal1/LTF/masked_genome/LTF.masked_genome.repeat_masked.fa'
window = 50000
nind = 40
theta = {}
outfile = '/mnt/gluster/home/sonal.singhal1/LTF/analysis/pop_gen/wattersons_theta.csv'
vcfbase = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.%s.filtered.coverage.repeatmasked.recoded_biallelicSNPs.vqsr.vcf.gz'
#vcfbase = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.%s.coverage.repeatmasked.filtered.recoded_biallelicSNPs.nomendel.vcf.gz'

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


for chr in chr_lengths:
	gzvcf = vcfbase % chr
	f = gzip.open(gzvcf, 'r')
	for l in f:
		if not re.search('^#', l):
			d = re.split('\t', l)
			# identifies the chunk into which the variant goes
			ix = int(int(d[1])/float(window))
			theta[chr][ix]['ss'] += 1
	f.close()

harmonic = 0
for i in range(1, nind):
	harmonic += 1/float(i)

out = open(outfile, 'w')
out.write('chr,index,start,end,seq_length,watterson_theta\n')
for chr in theta:
	for ix in sorted(theta[chr].keys()):
		if theta[chr][ix]['seq_length'] > 0:
			thetaw = theta[chr][ix]['ss'] / float(theta[chr][ix]['seq_length']) / harmonic
			out.write('%s,%s,%s,%s,%s,%.4f\n' % (chr, ix, theta[chr][ix]['start'], theta[chr][ix]['end'], theta[chr][ix]['seq_length'], thetaw))
		else:
			out.write('%s,%s,%s,%s,0,NA\n' % (chr, ix, theta[chr][ix]['start'], theta[chr][ix]['end']))
out.close()
