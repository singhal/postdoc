import re

chr_lengths = { 'chr10': 20806668, 'chr11': 21403021, 'chr12': 21576510, 'chr13': 16962381,
                'chr14': 16419078, 'chr15': 14428146, 'chr16': 9909, 'chr17': 11648728,
                'chr18': 11201131, 'chr19': 11587733, 'chr1A': 73657157, 'chr1B': 1083483,
                'chr1': 118548696, 'chr20': 15652063, 'chr21': 5979137, 'chr22': 3370227,
                'chr23': 6196912, 'chr24': 8021379, 'chr25': 1275379, 'chr26': 4907541,
                'chr27': 4618897, 'chr28': 4963201, 'chr2': 156412533, 'chr3': 112617285,
                'chr4A': 20704505, 'chr4': 69780378, 'chr5': 62374962, 'chr6': 36305782,
                'chr7': 39844632, 'chr8': 27993427, 'chr9': 27241186, 'chrLG2': 109741,
                'chrLG5': 16416, 'chrLGE22': 883365, 'chrZ': 72861351 }
repeat_file = '/mnt/gluster/home/sonal.singhal1/reference/taeGut1.repeatMaskerBlast.repeatLibrary20140131.out'

chr_repeats = {}
for chr in chr_lengths:
	chr_repeats[chr] = 0
f = open(repeat_file, 'r')
for l in f:
	d = re.split('\s+', l.rstrip())
	# not a header line
	if re.search('^\d', d[0]):
		chr = d[4]
		start = int(d[5])
		end = int(d[6])
		if chr in chr_repeats:
			chr_repeats[chr]  += abs(end - start)
f.close()

print 'chromosome\tchromosome_length\tsummed_repeat_length\tfraction_repeats'
for chr in sorted(chr_repeats.keys()):
	print '%s\t%s\t%s\t%.5f' % (chr, chr_lengths[chr], chr_repeats[chr], (chr_repeats[chr] / float(chr_lengths[chr])))
