import re
import os
import subprocess

genome = '/mnt/gluster/home/sonal.singhal1/ZF/masked_genome/ZF.masked_genome.repeat_masked.fa'
switchfile = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/breaks/badintervals.%s.hapi.csv'
out = '/mnt/gluster/home/sonal.singhal1/ZF/masked_genome/ZF.masked_genome.repeat_masked.switch_masked.fa'

chrs = { 'chr10': 20806668, 'chr11': 21403021, 'chr12': 21576510, 'chr13': 16962381,
                'chr14': 16419078, 'chr15': 14428146, 'chr16': 9909, 'chr17': 11648728,
                'chr18': 11201131, 'chr19': 11587733, 'chr1A': 73657157, 'chr1B': 1083483,
                'chr1': 118548696, 'chr20': 15652063, 'chr21': 5979137, 'chr22': 3370227,
                'chr23': 6196912, 'chr24': 8021379, 'chr25': 1275379, 'chr26': 4907541,
                'chr27': 4618897, 'chr28': 4963201, 'chr2': 156412533, 'chr3': 112617285,
                'chr4A': 20704505, 'chr4': 69780378, 'chr5': 62374962, 'chr6': 36305782,
                'chr7': 39844632, 'chr8': 27993427, 'chr9': 27241186, 'chrLG2': 109741,
                'chrLG5': 16416, 'chrLGE22': 883365, 'chrZ': 72861351}

switch = {}
for chr in chrs:
	f = open(switchfile % chr, 'r')
	header = f.next()
	for l in f:
		d = re.split(',', l.rstrip())
		if chr not in switch:
			switch[chr] = {}
		start = int(d[0]) + 1
		end = int(d[1])

		for pos in range(start, end):
			switch[chr][pos] = 1
	f.close()

def get_chromosome(genome, chr):
        outfile = genome + '_' + chr
        subprocess.call('~/bin/samtools-0.1.19/samtools faidx %s %s > %s' % (genome, chr, outfile), shell=True)
        out_f = open(outfile, 'r')
        chromosome = ''
        locus_name = out_f.next()
        for l in out_f:
                chromosome = chromosome + l.rstrip()
        out_f.close()
        os.remove(outfile)
        return list(chromosome)

out_f = open(out, 'w')
for chr in chrs:
	chr_as_list = get_chromosome(genome, chr)
	if chr in switch:
		for pos in switch[chr]:
			chr_as_list[pos - 1] = '8'
	out_f.write('>%s\n' % (chr))
	for i in xrange(0, len(chr_as_list), 60):
		out_f.write(''.join(chr_as_list[i:i+60]) + '\n')
out_f.close()
