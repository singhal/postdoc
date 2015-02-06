import re
import glob

chr_lengths = { 'chr10': 20806668, 'chr11': 21403021, 'chr12': 21576510, 'chr13': 16962381,
                'chr14': 16419078, 'chr15': 14428146, 'chr16': 9909, 'chr17': 11648728,
                'chr18': 11201131, 'chr19': 11587733, 'chr1A': 73657157, 'chr1B': 1083483,
                'chr1': 118548696, 'chr20': 15652063, 'chr21': 5979137, 'chr22': 3370227,
                'chr23': 6196912, 'chr24': 8021379, 'chr25': 1275379, 'chr26': 4907541,
                'chr27': 4618897, 'chr28': 4963201, 'chr2': 156412533, 'chr3': 112617285,
                'chr4A': 20704505, 'chr4': 69780378, 'chr5': 62374962, 'chr6': 36305782,
                'chr7': 39844632, 'chr8': 27993427, 'chr9': 27241186, 'chrLG2': 109741,
                'chrLG5': 16416, 'chrLGE22': 883365, 'chrZ': 72861351}


print 'chromosome,length_dropped,chr_length,percent_dropped'
files = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/breaks/bad*')
for file in files:
	chr = re.search('(chr[A-Z|0-9]+)', file).group(1)
	total = 0
	o = open(file, 'r')
	header = o.next()
	for l in o:
		d = re.split(',', l.rstrip())
		total += (int(d[1]) - int(d[0]) + 1)
	print '%s,%s,%s,%.3f' % (chr, total, chr_lengths[chr], float(total) / chr_lengths[chr])
