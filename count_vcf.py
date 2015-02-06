import re
import subprocess
import glob

chr_lengths = { 'chr10': 20806668, 'chr11': 21403021, 'chr12': 21576510, 'chr13': 16962381,
                'chr14': 16419078, 'chr15': 14428146, 'chr16': 9909, 'chr17': 11648728,
                'chr18': 11201131, 'chr19': 11587733, 'chr1A': 73657157, 'chr1B': 1083483,
                'chr1': 118548696, 'chr20': 15652063, 'chr21': 5979137, 'chr22': 3370227,
                'chr23': 6196912, 'chr24': 8021379, 'chr25': 1275379, 'chr26': 4907541,
                'chr27': 4618897, 'chr28': 4963201, 'chr2': 156412533, 'chr3': 112617285,
                'chr4A': 20704505, 'chr4': 69780378, 'chr5': 62374962, 'chr6': 36305782,
                'chr7': 39844632, 'chr8': 27993427, 'chr9': 27241186, 'chrLG2': 109741,
                'chrLG5': 16416, 'chrLGE22': 883365, 'chrZ': 72861351 }

files = glob.glob('/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/for_shapeit/*ltf*')

for file in files:
	chr = re.search('(chr[0-9|A-Z]+)', file).group(1)

	count = subprocess.Popen('zcat %s | tail -n 1' % file, shell=True, stdout=subprocess.PIPE)
	for l in count.stdout:
		d = re.split('\t', l)
		last_snp = int(d[1])
		div = last_snp / float(chr_lengths[chr])
	
		print '%s %s %s %.3f' % (chr, last_snp, chr_lengths[chr], div)
