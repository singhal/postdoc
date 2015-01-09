import re
from itertools import izip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--chr", help="chromosome for which to run analysis")
args = parser.parse_args()
chr = args.chr

dir = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/'
out = '%sbreaks/recombinationbreaks_across_lengths.%s.hapi.csv' % (dir, chr)

o = open(out, 'w')
o.write('chromosome,breaks,lengths,all_sites,het_sites\n')

hap_file = '%soutput_files/hapi.%s.csv' % (dir, chr)

# store haplotypes
haps = {}
f = open(hap_file, 'r')
for l in f:
	d = re.split(',', l.rstrip())
	hap = d[1:-1]
	if d[0] not in haps:
		haps[d[0]] = []
	haps[d[0]].append(hap)
f.close()

children = ['21', '22', '23']
parents = ['12', '10']

# identify heterozygous sites in child
for child in children:
	# paternal haplotype always comes first
	for parent, child_hap_wmiss in zip(parents, haps[child]):	
		match = []
		for ix, (a, b, c) in enumerate(izip(child_hap_wmiss, haps[parent][0], haps[parent][1])):
			if a != '0' and b != '0' and c != '0':
				if b != c:
					if a == b:
						match.append('0')
					elif a == c:
						match.append('1')
		
		matches = ''.join(match)
		# first do it for files as is
		switch = 0
                for ix, (pos1, pos2) in enumerate(zip(matches, matches[1:])):
                        if pos1 != pos2:
                                switch += 1    
                o.write('%s,%s,%s,%s,%s\n' % (chr, switch, '0', len(child_hap_wmiss), len(match)))

		for bp in range(1,10):
			matches = re.sub('01{%s}0' % bp, '00', matches)
			matches = re.sub('10{%s}1' % bp, '11', matches)
						
			switch = 0
			for ix, (pos1, pos2) in enumerate(zip(matches, matches[1:])):
				if pos1 != pos2:
					switch += 1	
			o.write('%s,%s,%s,%s,%s\n' % (chr, switch, bp, len(child_hap_wmiss), len(match)))
	
		for bp in [100, 500, 1000, 5000, 10000, 50000]:
			matches = re.sub('01{1,%s}0' % bp, '00', matches)
                        matches = re.sub('10{1,%s}1' % bp, '11', matches)
                                                
                        switch = 0
                        for ix, (pos1, pos2) in enumerate(zip(matches, matches[1:])):
                                if pos1 != pos2:
                                        switch += 1
                        o.write('%s,%s,%s,%s,%s\n' % (chr, switch, bp, len(child_hap_wmiss), len(match)))
o.close()
