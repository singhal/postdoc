import re
import subprocess
import numpy as np

chr_lengths = { 'chr10': 20806668, 'chr11': 21403021, 'chr12': 21576510, 'chr13': 16962381,
                'chr14': 16419078, 'chr15': 14428146, 'chr16': 9909, 'chr17': 11648728,
                'chr18': 11201131, 'chr19': 11587733, 'chr1A': 73657157, 'chr1B': 1083483,
                'chr1': 118548696, 'chr20': 15652063, 'chr21': 5979137, 'chr22': 3370227,
                'chr23': 6196912, 'chr24': 8021379, 'chr25': 1275379, 'chr26': 4907541,
                'chr27': 4618897, 'chr28': 4963201, 'chr2': 156412533, 'chr3': 112617285,
                'chr4A': 20704505, 'chr4': 69780378, 'chr5': 62374962, 'chr6': 36305782,
                'chr7': 39844632, 'chr8': 27993427, 'chr9': 27241186, 'chrLG2': 109741,
                'chrLG5': 16416, 'chrLGE22': 883365, 'chrZ': 72861351, 'chrZ_random': 500000}
window = 50000
genome = '/mnt/gluster/home/sonal.singhal1/reference/Zfinch.fa'
out = '/mnt/gluster/home/sonal.singhal1/reference/reference.cpg_cg.csv'
o = open(out, 'w')
o.write('chr,ix,start,end,num_sites,cpg,gc\n')

for chr, length in chr_lengths.items():
        for ix, start in enumerate(range(1, length, window)):
                end = start + window - 1
                if end > length:
                        end = length
               
		seq_call = subprocess.Popen('samtools faidx %s %s:%s-%s' % (genome, chr, start, end), shell=True, stdout=subprocess.PIPE)
             	seq = ''
                for l in seq_call.stdout:
                        if not re.match('>', l):
				seq += l.rstrip().upper()
		
		cpg = 'NA'
		tuples = [seq[i:i+2] for i in range(0, len(seq), 2)] + [seq[i:i+2] for i in range(1, len(seq), 2)]
		tuples = [x for x in tuples if not re.search('N', x)]
		if len(tuples) > 0:
			cpg = (tuples.count('CG') * 2) / float(len(tuples))

		seq = list(seq)
		counts = {}
		for nuc in ['A', 'T', 'C', 'G']:
			counts[nuc] = seq.count(nuc)
		gc = 'NA'
		if np.sum(counts.values()) > 0:
			gc = (counts['C'] + counts['G']) / float( np.sum(counts.values()) )

		o.write('%s,%s,%s,%s,%s,%s,%s\n' % (chr, ix, start, end, np.sum(counts.values()), cpg, gc))
o.close()
