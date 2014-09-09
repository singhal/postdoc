import pandas as pd
import numpy as np
import re
from itertools import izip
import subprocess
import random

def get_data(file, seq_size, rho_file, num_sim):
	d = pd.read_csv(file, sep=" ", skiprows=3, header=None,
                        names=['left_snp', 'right_snp', 'meanrho', 'p0.025', 'p0.5', 'p0.975'])
	
	# chunk the sequence into the desired simulation length
	bins = range(d.left_snp.min(), d.left_snp.max(), seq_size)
	groups = d.groupby(np.digitize(d.left_snp, bins))

	# create a hash of the median recombination rate for each simulation length sequence 
	rec_rates = {}
	for ix, group in groups:
	        rec_rates[bins[ix - 1]] = group.meanrho.median()
	ordered_chunks = sorted(rec_rates, key=rec_rates.get)

	# figure out which of the chunks will be accessed for the simulations
	# do this all because want a spread of recombination rates
	# should rewrite this using the percentile function
	sim = [0] + [int(len(ordered_chunks) * float(i)/ num_sim) for i in  range(1, num_sim + 1) ]
	sim[-1] = sim[-1] - 1	

	rho_f = open(rho_file, 'w')
	for ix, sim_num in enumerate(sim):
		rho_f.write('%s %s %s\n' % (ix, ordered_chunks[sim_num], rec_rates[ordered_chunks[sim_num]]))
	rho_f.close()

	return d, sim, ordered_chunks, rec_rates

# this is the worst function. i really need to fix it
def create_simulation_files(ix, i, d, sim, ordered_chunks, rec_rates, out_dir, seq_size, num_rep, chunk_size, diff, theta, nsam, eq_freq, mut_rates):
	chunk_start = ordered_chunks[i]
	chunk_end = chunk_start + seq_size

	tmp = d[ d.left_snp <= chunk_end ]
	tmp = tmp[ tmp.right_snp >= chunk_start ]

	chunk_rec = {}
	for rho, left, right in izip(d.meanrho, d.left_snp, d.right_snp):
		for bp in range(left, right):
			chunk_rec[ bp ] = rho
	
	median_rate = rec_rates[ordered_chunks[i]]
	spots = {}
	chunk_kb = {}
	for kb1_start in range(int(chunk_start), int(chunk_end), chunk_size):
		kb1_end = kb1_start + chunk_size
		sum = 0
		for bp in range(kb1_start, kb1_end):
			sum += chunk_rec[bp]
		kb_avg = sum / float(chunk_size)
		if kb_avg / median_rate > diff or kb_avg / median_rate < 1/float(diff):
			spots[ kb1_start ] = kb_avg
	
	hotspot_file = '%shotspot%s.txt' % (out_dir, ix)
	hotspot_f = open(hotspot_file, 'w')
	for spot in sorted(spots.keys()):
		seq_frac_start = (spot - chunk_start) / float(seq_size)
		seq_frac_end = (spot + chunk_size - 1 - chunk_start) / float(seq_size)
		hotspot_f.write('%.4f\t%.4f\t%.5f\n' % (seq_frac_start, seq_frac_end, spots[spot] / median_rate))
	hotspot_f.close()	

	for j in range(num_rep):
		haplo_file = '%shaplo%s_%s.fa' % (out_dir, ix, j)
		anc_file = '%sancallele%s_%s.txt' % (out_dir, ix, j)
		haplo_f = open(haplo_file, 'w')
		anc_f = open(anc_file, 'w')

		macs = subprocess.Popen('~/bin/macs/macs %s %s -t %s -r %s -R %s | ~/bin/macs/msformatter' % (nsam, seq_size, theta, median_rate, hotspot_file), shell=True, stdout=subprocess.PIPE)

		positions = []
		haplo = []
		for l in macs.stdout:
			if re.match('^positions:', l):
				positions = [float(match) for match in re.findall('([\d\.e-]+)', l)]
			haplo.append(l.rstrip())			

		haplo = haplo[(len(haplo) - nsam):len(haplo)]

		bases = []
		for base, freq in eq_freq.items():
			bases += [base] * int(freq * 1000)
		seq = []
		for x in range(int(seq_size)):
			seq.append(random.choice(bases))

		mutations = {}
	
		positions = [int(round(pos * seq_size)) - 1 for pos in positions]
	
		for pos_ix, pos in enumerate(positions):
			if pos < 0:
				positions[pos_ix] = 0
			if pos > (seq_size -1):
				positions[pos_ix] = seq_size - 1
		for pos in positions:
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

		for ind, ind_hap in enumerate(haplo):
			tmp_seq = seq[:]
			for hap_ix, bp in enumerate(list(ind_hap)):
				if bp == '1':
					tmp_seq[positions[hap_ix]] = mutations[positions[hap_ix]][1]
			haplo_f.write('>haplo%s\n' % ind)
			for seq_i in xrange(0, len(tmp_seq), 60):
				haplo_f.write('%s\n' % ''.join(tmp_seq[seq_i:seq_i+60]))
		haplo_f.close()


def main():
	out_dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/simulations/'
	file = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/drafts_Aug30/chr1A.bpen10.txt'
	# sim MB
	seq_size = 1000000
	# num MB to simulate
	num_sim = 11
	# num replicates
	num_rep = 3
	# smoothing parameter
	chunk_size = 1000
	# hotspot / coldspot difference, > 1
	diff = 3
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
	
	rho_file = '%smean_rho.txt' % (out_dir)
	d, sim, ordered_chunks, rec_rates = get_data(file, seq_size, rho_file, num_sim)
	for ix, i in enumerate(sim):
		create_simulation_files(ix, i, d, sim, ordered_chunks, rec_rates, out_dir, seq_size, num_rep, chunk_size, diff, theta, nsam, eq_freq, mut_rates)

if __name__ == "__main__":
    main()

