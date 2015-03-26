import numpy as np
import pandas as pd
import re
import subprocess
from itertools import izip
import pickle
import copy


def get_hotspots(hotspot_file):
	hot = pd.read_csv(hotspot_file)
        hot = hot[hot.zlk >= 10]

	hotspots = {}

	for chr, start, length in izip(hot.chr, hot.zstart, hot.zlength):
		if chr not in hotspots:
			hotspots[chr] = {}
		midpoint = start + int(length / 2.0)
		hotspots[chr][midpoint] = {'GC': None, 'CpG': None, 'match': None, 'length': length}

	return hotspots


def get_coldspots(coldspot_file):
	d = pd.read_csv(coldspot_file)
	coldspots = {}
	for chr, start, end in izip(d.chr, d.spot_start, d.spot_end):
		if chr not in coldspots:
			coldspots[chr] = {}
		midpoint = int((start + end)/2.0)
		coldspots[chr][midpoint] = {'GC': None, 'CpG': None, 'match': None, 'length': (end - start)}
	
	return coldspots


def get_cg_data(genome, spots):
	for chr in spots:
		for center in spots[chr]:
			start = center - 500
			end = center + 500

			seq = ''
			out = subprocess.Popen('samtools faidx %s %s:%s-%s' % (genome, chr, start, end), shell=True, stdout=subprocess.PIPE)
                	for l in out.stdout:
                        	if not re.match('>', l):
                                	seq += l.rstrip().upper()
	
			tuples = [seq[i:i+2] for i in range(0, len(seq), 2)] + [seq[i:i+2] for i in range(1, len(seq), 2)]
			tuples = [x for x in tuples if not re.search('N', x)]
			if len(tuples) > 0:
				spots[chr][center]['CpG'] = tuples.count('CG') / float(len(tuples))

			seq = list(seq)	
			counts = {}
			for nuc in ['A', 'T', 'C', 'G']:
				counts[nuc] = seq.count(nuc)
			if np.sum(counts.values()) > 0:
				spots[chr][center]['GC'] = (counts['C'] + counts['G']) / float( np.sum(counts.values()) )
	
	return spots			


def match_hotspots(hot, cold, cpg_limit, gc_limit):
	for chr in hot:
		for hot_center in hot[chr]:
			gc_hot = hot[chr][hot_center]['GC']
			cpg_hot = hot[chr][hot_center]['CpG']

			if gc_hot and cpg_hot:
				possible_matches = []
				for center in cold[chr]:
					# don't want to consider spots that are already matched
					if not cold[chr][center]['match']:
						gc_cold = cold[chr][center]['GC']
                        			cpg_cold = cold[chr][center]['CpG']
				
						if gc_cold and cpg_cold:
							if abs(gc_cold - gc_hot) <= gc_limit and abs(cpg_cold - cpg_hot) <= gc_hot:
								possible_matches.append(center)

				dist = 1e9
				match_cold = None
				for match in possible_matches:
					if abs(match - hot_center) < dist:
						dist = abs(match - hot_center)
						match_cold = match
				if match_cold:
					hot[chr][hot_center]['match'] = match_cold
					cold[chr][match_cold]['match'] = hot_center
				
	return hot, cold


def rev_comp(word):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
	bases = list(word) 
    	bases = [complement[base] for base in bases] 

	return ''.join(bases[::-1])


def get_kmers(kmer_min, kmer_max):
	kmers = {}

	for k in range(kmer_min - 1, kmer_max):
		bases = ['A','C','G','T']
		words = copy.copy(bases)

		for i in range(k):
			newwords = []
			for w in words:
				for b in bases:
					newwords.append(w + b)
			words = []
			words = newwords

		for word in words:
		     if word not in kmers and rev_comp(word) not in kmers:
			kmers[ word ] = {'hot': 0, 'cold': 0}

	return kmers


def test_seq(genome, spots, type, kmers):
	
	for chr in spots:
		for center in spots[chr]:
			if spots[chr][center]['match']:
				start = center - 1500
                        	end = center + 1500

                        	seq = ''
                        	out = subprocess.Popen('samtools faidx %s %s:%s-%s' % (genome, chr, start, end), shell=True, stdout=subprocess.PIPE)
                        	for l in out.stdout:
                        	        if not re.match('>', l):
                                	        seq += l.rstrip().upper()

				for kmer in kmers:
					kmer1 = kmer
					kmer2 = rev_comp(kmer)

					match = False
					if seq.find(kmer1) > -1:
						match = True
					if seq.find(kmer2) > - 1:
						match = True
			
					if match:
						kmers[kmer][type] += 1

	return kmers


def main():
	genome = '/mnt/gluster/home/sonal.singhal1/reference/ancestral_genome.fa'
	out_dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/motifs/'
	cpg_limit = 0.001
	gc_limit = 0.005
	kmer_min = 3
	kmer_max = 10

	# first match hotspots and coldspots
	# get the hotspots
	hotspot_file = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/spot2kb_flank40kb.seqldhot_validate_hotspots.csv'
	hot = get_hotspots(hotspot_file)
	# get the coldspots
	coldspot_file = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/ZF.ldhelmet_unvalidated_coldspots.csv'
	cold = get_coldspots(coldspot_file)
	# then get cpg / CG data
	hot = get_cg_data(genome, hot)
	cold = get_cg_data(genome, cold)

	pickle.dump(hot, open('%shotspots.pickle' % out_dir, 'w'))
	pickle.dump(cold, open('%scoldspots.pickle' % out_dir, 'w'))
	
	# hot = pickle.load(open('hotspots.pickle', 'r'))
	# cold = pickle.load(open('coldspots.pickle', 'r'))
	
	# hot, cold = match_hotspots(hot, cold, cpg_limit, gc_limit)	
	# pickle.dump(hot, open('hotspots.pickle', 'w'))
        # pickle.dump(cold, open('coldspots.pickle', 'w'))

	# then generate kmer hash -- make it hash[ kmer ] = {'cold': #, 'hot': #}
	# kmers = get_kmers(kmer_min, kmer_max)
	
	# then test each hotspot and coldspot
	# kmers = test_seq(genome, hot, 'hot', kmers)
	# kmers = test_seq(genome, cold, 'cold', kmers)

	# then print out results
	# for kmer in kmers:
	#	print '%s,%s,%s' % (kmer, kmers[kmer]['hot'], kmers[kmer]['cold'])

if __name__ == "__main__":
    main()
