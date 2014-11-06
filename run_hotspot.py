import re
import sys
import os

chrs = { 'chr10': 20806668, 'chr11': 21403021, 'chr12': 21576510, 'chr13': 16962381,
                'chr14': 16419078, 'chr15': 14428146, 'chr16': 9909, 'chr17': 11648728,
                'chr18': 11201131, 'chr19': 11587733, 'chr1A': 73657157, 'chr1B': 1083483,
                'chr1': 118548696, 'chr20': 15652063, 'chr21': 5979137, 'chr22': 3370227,
                'chr23': 6196912, 'chr24': 8021379, 'chr25': 1275379, 'chr26': 4907541,
                'chr27': 4618897, 'chr28': 4963201, 'chr2': 156412533, 'chr3': 112617285,
                'chr4A': 20704505, 'chr4': 69780378, 'chr5': 62374962, 'chr6': 36305782,
                'chr7': 39844632, 'chr8': 27993427, 'chr9': 27241186, 'chrLG2': 109741,
                'chrLG5': 16416, 'chrLGE22': 883365, 'chrZ': 72861351}

tracker = 1
for species in ['ZF', 'LTF']:
	for bp in [500, 1000, 2000]:
		for flank in [20000, 40000]:
	#for bp in [1000]:
	#	for flank in [40000]:
			for chr in chrs:
				file = '/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/putative_hotspots/%s.putativehotspots.block%s_flank%s.out' % (species, chr, bp, flank)
				if not os.path.isfile(file):
					if not (chr == 'chr16' and species == 'LTF'):
						call = "python ~/scripts/find_hotspots.py --chr \'%s\' --block %s --flank %s --sp \'%s\'" % (chr, bp, flank, species)
						print 'echo \"%s\" | qsub -l h_vmem=5g -cwd -V -j y -N \'find%s\'' % (call, tracker)
						#print file
						tracker += 1
