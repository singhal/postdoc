import glob
import re
import numpy as np

files = glob.glob('/mnt/gluster/home/sonal.singhal1/gene_trees/individual_trees/phyml_files/*stats*')

base_freq = {'A':[], 'T': [], 'G': [] , 'C': [] }
rel_rate = {'AC': [], 'AG': [], 'AT': [], 'CG': [], 'CT': [], 'GT': []}

for file in files:
	f = open(file, 'r')
	for l in f:
		if re.search('f\(\S\)=', l):
			match = re.search('f\((\S)\)=\s+([\d\.]+)', l)
			base_freq[match.group(1)].append(float(match.group(2)))
		if re.search('\S\s+<->\s+\S', l):
			match = re.search('(\S)\s+<->\s+(\S)\s+([\d\.]+)', l)
			rel_rate['%s%s' % (match.group(1), match.group(2))].append(float(match.group(3)))
	f.close()

for base, freqlist in base_freq.items():
	print '%s %.3f' % (base, np.mean(freqlist))
for transtype, ratelist in rel_rate.items():
	print '%s %.3f' % (transtype, np.mean(ratelist))
		
