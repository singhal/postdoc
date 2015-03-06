import re
import gzip
import glob
import subprocess
import argparse
import random

parser = argparse.ArgumentParser()
parser.add_argument("--sp", help="species for which to run analysis")
args = parser.parse_args()
sp = args.sp

chr_lengths = { 'chr10': 20806668, 'chr11': 21403021, 'chr12': 21576510, 'chr13': 16962381,
                'chr14': 16419078, 'chr15': 14428146, 'chr16': 9909, 'chr17': 11648728,
                'chr18': 11201131, 'chr19': 11587733, 'chr1A': 73657157, 'chr1B': 1083483,
                'chr1': 118548696, 'chr20': 15652063, 'chr21': 5979137, 'chr22': 3370227,
                'chr23': 6196912, 'chr24': 8021379, 'chr25': 1275379, 'chr26': 4907541,
                'chr27': 4618897, 'chr28': 4963201, 'chr2': 156412533, 'chr3': 112617285,
                'chr4A': 20704505, 'chr4': 69780378, 'chr5': 62374962, 'chr6': 36305782,
                'chr7': 39844632, 'chr8': 27993427, 'chr9': 27241186, 'chrLG2': 109741,
                'chrLG5': 16416, 'chrLGE22': 883365, 'chrZ': 72861351}

if sp == 'LTF':
        genome = '/mnt/gluster/home/sonal.singhal1/LTF/masked_genome/LTF.masked_genome.repeat_masked.fa'
if sp == 'ZF':
        genome = '/mnt/gluster/home/sonal.singhal1/ZF/masked_genome/ZF.masked_genome.repeat_masked.switch_masked.fa'

window = 50000

if sp == 'ZF':
        vcfs = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/*phased*')
if sp == 'LTF':
        vcfs = glob.glob('/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/*phased*')

pi = {}
outfile = '/mnt/gluster/home/sonal.singhal1/%s/analysis/pop_gen/pi.csv' % sp
o = open(outfile, 'w')
o.write('chr,index,start,end,seq_length,pi\n')

for chr, length in chr_lengths.items():
        pi[chr] = {}
        for ix, start in enumerate(range(1, length, window)):
                end = start + window - 1
                if end > length:
                        end = length
                
                seq = subprocess.Popen('samtools faidx %s %s:%s-%s' % (genome, chr, start, end), shell=True, stdout=subprocess.PIPE)
                n_count = 0
                for l in seq.stdout:
                        if not re.match('>', l):
                                n_count += len(filter(lambda x: int(x) < 4, list(l.rstrip())))

                pi[chr][ix] = {'start': start, 'end': end, 'num_diff': 0, 'seq_length': n_count}


for gzvcf in vcfs:
        f = gzip.open(gzvcf, 'r')
        chr = re.search('(chr[A-Z|0-9]+)', gzvcf).group(1)
	haps = {}
	nind = None
        for l in f:
                if not re.search('^#', l):
                        d = re.split('\t', l.rstrip())
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
					if re.search('/', geno):
						geno = re.split('/', geno)
						random.shuffle(geno)
						genos += geno
					else:
                                        	genos += re.split('\|', geno)
				haps[int(d[1])] = genos
				nind = len(genos)

	sites = sorted(haps.keys())	

	for ix in pi[chr]:
		start = pi[chr][ix]['start']
		end = pi[chr][ix]['end']

		tmp_sites = filter(lambda x: x >= start, sites)
		tmp_sites = filter(lambda x: x<= end, tmp_sites)

		tmp_haps = []
		for i in range(nind):
			tmp_hap = ''
			for site in tmp_sites:
				tmp_hap += haps[site][i]
			tmp_haps.append(tmp_hap)
		
		pi_sum = 0
		for i in range(nind):
			for j in range(i):
				# does this handle missing data correctly? no.
				diff = [x == y for (x, y) in zip(tmp_haps[i], tmp_haps[j])].count(False)
				if  pi[chr][ix]['seq_length'] > 0:
					pi_sum += (1/float(nind)) * (1/float(nind)) * (diff/float(pi[chr][ix]['seq_length']))
		if pi[chr][ix]['seq_length'] > 0:
			pi_sum = pi_sum * 2
			o.write('%s,%s,%s,%s,%s,%.4f\n' % (chr, ix, start, end, pi[chr][ix]['seq_length'], pi_sum))
		else:
			o.write('%s,%s,%s,%s,0,NA\n' % (chr, ix, start, end))
	
