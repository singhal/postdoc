import re
import os
import subprocess
import gzip

repeat = '/mnt/gluster/home/sonal.singhal1/reference/taeGut1.repeatMaskerBlast.repeatLibrary20140131.out'
vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/gatk.ug.all_zf.%s.coverage.filtered.vqsr.vcf.gz' 

chrs = { 'chr10': 20806668, 'chr11': 21403021, 'chr12': 21576510, 'chr13': 16962381,
                'chr14': 16419078, 'chr15': 14428146, 'chr16': 9909, 'chr17': 11648728,
                'chr18': 11201131, 'chr19': 11587733, 'chr1A': 73657157, 'chr1B': 1083483,
                'chr1': 118548696, 'chr20': 15652063, 'chr21': 5979137, 'chr22': 3370227,
                'chr23': 6196912, 'chr24': 8021379, 'chr25': 1275379, 'chr26': 4907541,
                'chr27': 4618897, 'chr28': 4963201, 'chr2': 156412533, 'chr3': 112617285,
                'chr4A': 20704505, 'chr4': 69780378, 'chr5': 62374962, 'chr6': 36305782,
                'chr7': 39844632, 'chr8': 27993427, 'chr9': 27241186, 'chrLG2': 109741,
                'chrLG5': 16416, 'chrLGE22': 883365, 'chrZ': 72861351}

repeats = {}
f = open(repeat, 'r')
for l in f:
	if re.search('^\s+\d+', l):
		d = re.split('\s+', l.rstrip())
		if d[5] in chrs:
			if d[5] not in repeats:
				repeats[d[5]] = {}
			start = int(d[6])
			end = int(d[7]) + 1

			for pos in range(start, end):
				repeats[d[5]][str(pos)] = 1
f.close()

for chr in chrs:
	vcffile = vcf % chr
	out = vcffile.replace('coverage.filtered', 'coverage.filtered.repeatmasked')
	
	f = gzip.open(vcffile, 'r')
	o = gzip.open(out, 'w')
	
	for l in f:
		l = l.rstrip()
		if re.search('^#', l):
			o.write(l + '\n')
		else:
			if re.search('PASS', l):
				d = re.split('\t', l)
				if d[0] in repeats:
					if d[1] not in repeats[d[0]]:
						o.write(l + '\n')
				else:
					o.write(l + '\n')
	f.close()
	o.close()
