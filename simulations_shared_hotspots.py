import numpy as np
import re
import sys
from itertools import izip
import subprocess
import random
import os
import argparse

def simulate(out_dir, seq_size, theta, nsam, eq_freq, mut_rates, rho, diffs, hotspot_lengths, ix):
	
	hotspot_file = '%shotspot_rho%s_diff%s_%s.txt' % (out_dir, rho, diffs[0], ix)
	hotspot_f = open(hotspot_file, 'w')
	num_hotspots = len(diffs) * len(hotspot_lengths)
	starts = np.linspace(0+50e3,seq_size-50e3,num_hotspots)
	for i, diff in enumerate(diffs):
		for j, length in enumerate(hotspot_lengths):
			spot_start = starts[j + i * len(hotspot_lengths)]
			spot_end = spot_start + length
			hotspot_f.write('%.4f\t%.4f\t%.5f\n' % (spot_start / float(seq_size), spot_end / float(seq_size), diff))
	hotspot_f.close()	

	haplo_file = '%shaplo_rho%s_diff%s_%s.fa' % (out_dir, rho, diffs[0], ix)
	anc_file = '%sancallele_rho%s_diff%s_%s.txt' % (out_dir, rho, diffs[0], ix)
	haplo_f = open(haplo_file, 'w')
	anc_f = open(anc_file, 'w')

	macs = subprocess.Popen('~/bin/macs/macs %s %s -t %s -r %s -R %s | ~/bin/macs/msformatter' % (nsam, seq_size, theta, rho, hotspot_file), shell=True, stdout=subprocess.PIPE)

	positions = []
	haplo = []
	for l in macs.stdout:
		if re.match('^positions:', l):
			positions = [float(match) for match in re.findall('([\d\.e-]+)', l)]
		haplo.append(l.rstrip())			

	haplo = haplo[(len(haplo) - nsam):len(haplo)]
	for ix, hap in enumerate(haplo):
		haplo[ix] = list(hap)
	# these don't have singletons
	singletons = []
	positions_pruned = []
	haplo_pruned = []
	for ix, pos in enumerate(positions):
		allele_count = 0
		for hap in haplo:
			if hap[ix] == '1':
				allele_count += 1
		if allele_count <= 1:
			singletons.append(ix)

	for ix, pos in enumerate(positions):
		if ix not in singletons:
			pos = int(round(pos * seq_size)) - 1 
			if pos < 0:
	                        pos = 0
        	        if pos > (seq_size -1):
                	        pos = seq_size - 1
			positions_pruned.append(pos)
	
	for ix1, hap in enumerate(haplo):
		haplo_pruned.append([])
		for ix2, pos in enumerate(hap):
			if ix2 not in singletons:
				haplo_pruned[ix1].append(pos)

	del positions
	del haplo
	del singletons
	
	bases = []
	for base, freq in eq_freq.items():
		bases += [base] * int(freq * 1000)
	seq = []
	for x in range(int(seq_size)):
		seq.append(random.choice(bases))

	mutations = {}
	
	for pos in positions_pruned:
		anc = seq[pos]
		ran_num = random.random()
		for base, prob in mut_rates[anc].items():
			ran_num = ran_num - prob
			if ran_num < 0:
				mut = base
				nuc = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
				prob = [0.02, 0.02, 0.02, 0.02]
				prob[nuc[mut]] = 0.94
				anc_f.write('%s %s\n' % (pos, ' '.join('%.2f' % i for i in prob)))
				mutations[pos] = [anc, mut]
				break
	anc_f.close()

	for ind, ind_hap in enumerate(haplo_pruned):
		tmp_seq = seq[:]
		for hap_ix, bp in enumerate(ind_hap):
			if bp == '1':
				tmp_seq[positions_pruned[hap_ix]] = mutations[positions_pruned[hap_ix]][1]
		haplo_f.write('>haplo%s\n' % ind)
		for seq_i in xrange(0, len(tmp_seq), 60):
			haplo_f.write('%s\n' % ''.join(tmp_seq[seq_i:seq_i+60]))
	haplo_f.close()


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--sp", help="species for which to run analysis")
	args = parser.parse_args()
	sp = args.sp

	if sp == 'ZF':
		out_dir = '/mnt/gluster/home/sonal.singhal1/simulations/shared/ZF/'
		theta = 0.0131
		nsam = 38
		#rhos = [0.8]
		rhos = [0.001, 0.002, 0.01, 0.1, 0.5, 0.8]
	if sp == 'LTF':
		out_dir = '/mnt/gluster/home/sonal.singhal1/simulations/shared/LTF/'
                theta = 0.0073
                nsam = 40
		#rhos = [0.0005, 0.001, 0.005, 0.4]
		rhos = [0.0005, 0.001, 0.005, 0.05, 0.25, 0.4]
	# sim MB
	seq_size = 1000000
	# num replicates to simulate
	num_sim = 10
	# hotspot / coldspot difference, > 1
	diffs_array = [[5]*12, [10]*12, [15]*12, [20]*12, [40]*12, [60]*12, [80]*12, [100]*12]
	# hotspot length
	hotspot_lengths = [2000]
	# A, C, T, G
	eq_freq = {'A': 0.303,'C': 0.197, 'G': 0.305, 'T': 0.195}
	# mutation rates, modified from mutation matrix
	mut_rates = {'A': {'C': 0.191, 'G': 0.591, 'T': 0.218},
                     'C': {'A': 0.206, 'G': 0.135, 'T': 0.659},
                     'G': {'A': 0.659, 'C': 0.135, 'T': 0.206},
                     'T': {'A': 0.215, 'C': 0.600, 'G': 0.185}}
	
	for rho in rhos:
		for ix in range(num_sim):
			for diffs in diffs_array:
				haplo_file = '%shaplo_rho%s_diff%s_%s.fa' % (out_dir, rho, diffs[0], ix)
				haplo_file1 = '%shaplo/haplo_rho%s_diff%s_%s.fa' % (out_dir, rho, diffs[0], ix)
				if not os.path.isfile(haplo_file) and not os.path.isfile(haplo_file1):
					print haplo_file
					#simulate(out_dir, seq_size, theta, nsam, eq_freq, mut_rates, rho, diffs, hotspot_lengths, ix)

if __name__ == "__main__":
    main()
