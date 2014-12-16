import numpy as np
import re
import sys
from itertools import izip
import subprocess
import random
import os


def switch_error(point, hap1, hap2):
	hap1new = ''.join(hap1[0:point]) + ''.join(hap2[point:])
	hap2new = ''.join(hap2[0:point]) + ''.join(hap1[point:])

	return list(hap1new), list(hap2new)


def simulate(out_dir, seq_size, theta, nsam, eq_freq, mut_rates, rho, switch, ix):
	
	haplo_file = '%shaplo_rho%s_switch%s_%s.fa' % (out_dir, rho, switch, ix)
	anc_file = '%sancallele_rho%s_switch%s_%s.txt' % (out_dir, rho, switch, ix)
	haplo_f = open(haplo_file, 'w')
	anc_f = open(anc_file, 'w')

	macs = subprocess.Popen('~/bin/macs/macs %s %s -t %s -r %s | ~/bin/macs/msformatter' % (nsam, seq_size, theta, rho), shell=True, stdout=subprocess.PIPE)

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

	random.shuffle(haplo_pruned)
	switch_hap = []
	# introduce switch error rate	
	for i in range(0, nsam, 2):
		hap1 = list(haplo_pruned[i])
		hap2 = list(haplo_pruned[i + 1])

		num_hets = []
		for pos, (i, j) in enumerate(zip(hap1, hap2)):
			if i != j:
				num_hets.append(pos)
		switch_points = sorted(random.sample(num_hets, int(len(num_hets) * switch)))

		for switch_point in switch_points:
			hap1, hap2 = switch_error(switch_point, hap1, hap2)
		switch_hap.append(hap1)
		switch_hap.append(hap2)
	del haplo_pruned
	
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

	for ind, ind_hap in enumerate(switch_hap):
		tmp_seq = seq[:]
		for hap_ix, bp in enumerate(ind_hap):
			if bp == '1':
				tmp_seq[positions_pruned[hap_ix]] = mutations[positions_pruned[hap_ix]][1]
		haplo_f.write('>haplo%s\n' % ind)
		for seq_i in xrange(0, len(tmp_seq), 60):
			haplo_f.write('%s\n' % ''.join(tmp_seq[seq_i:seq_i+60]))
	haplo_f.close()


def main():
	out_dir = '/mnt/gluster/home/sonal.singhal1/simulations/FDR/'
	# sim MB
	seq_size = 1000000
	# num replicates to simulate
	num_sim = 4
	# mean rho values
	rhos = [0.001, 0.01, 0.1, 0.5]
	# switch error rate
	switches = [0.0, 0.01, 0.02, 0.04, 0.08, 0.16]
	# theta (per bp!)
	theta = 0.0135
	# number of haplotypes to sample
	nsam = 38
	# A, C, T, G
	eq_freq = {'A': 0.303,'C': 0.197, 'G': 0.305, 'T': 0.195}
	# mutation rates, modified from mutation matrix
	mut_rates = {'A': {'C': 0.191, 'G': 0.591, 'T': 0.218},
                     'C': {'A': 0.206, 'G': 0.135, 'T': 0.659},
                     'G': {'A': 0.659, 'C': 0.135, 'T': 0.206},
                     'T': {'A': 0.215, 'C': 0.600, 'G': 0.185}}
	

	for rho in rhos:
		for switch in switches:
			for ix in range(num_sim):
				haplo_file = '%shaplotypes/haplo_rho%s_switch%s_%s.fa' % (out_dir, rho, switch, ix)
				if not os.path.isfile(haplo_file):
					simulate(out_dir, seq_size, theta, nsam, eq_freq, mut_rates, rho, switch, ix)
	

if __name__ == "__main__":
    main()

