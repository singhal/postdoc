import numpy as np
import pandas as pd
import re
import subprocess
from itertools import izip
import pickle
import copy


def get_hotspots(hotspot_file, sp):
	hot = pd.read_csv(hotspot_file)

        hot = hot[hot.llk >= 10]
	hot = hot[hot.zlk >= 10]

	hotspots = {}

	for chr, start, length in izip(hot.chr, hot.lstart, hot.llength):
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
			start1 = center - 5500
			end1 = center - 500
			
			start2 = center + 500
			end2 = center + 5500

			seq = ''
			out = subprocess.Popen('samtools faidx %s %s:%s-%s' % (genome, chr, start1, end1), shell=True, stdout=subprocess.PIPE)
                	for l in out.stdout:
                        	if not re.match('>', l):
                                	seq += l.rstrip().upper()

			out = subprocess.Popen('samtools faidx %s %s:%s-%s' % (genome, chr, start2, end2), shell=True, stdout=subprocess.PIPE)
                        for l in out.stdout:
                                if not re.match('>', l):
                                        seq += l.rstrip().upper()
	
			tuples = [seq[i:i+2] for i in range(0, len(seq), 2)] + [seq[i:i+2] for i in range(1, len(seq), 2)]
			tuples = [x for x in tuples if not re.search('N', x)]
			if len(tuples) > 0:
				spots[chr][center]['CpG'] = (tuples.count('CG') * 2) / float(len(tuples))

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
							if abs(gc_cold - gc_hot) <= gc_limit and abs(cpg_cold - cpg_hot) <= cpg_limit:
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


def print_hotspots(out_dir, hot, cold):
	out = '%smatched_hotspots.01GC.001CpG.all_finches.csv' % out_dir
	o = open(out, 'w')

	o.write('chr,hot_mid,hot_GC,hot_CpG,cold_mid,cold_GC,cold_CpG\n')
	for chr in hot:
		for hot_mid in hot[chr]:
			if hot[chr][hot_mid]['match']:
				cold_mid = hot[chr][hot_mid]['match']
				o.write('%s,%s,%s,%s,%s,%s,%s\n' % (chr, hot_mid, hot[chr][hot_mid]['GC'], hot[chr][hot_mid]['CpG'], 
									cold_mid, cold[chr][cold_mid]['GC'], cold[chr][cold_mid]['CpG']))
	o.close()


def main():
	genome = '/mnt/gluster/home/sonal.singhal1/reference/ancestral_genome.all_finches.fa'
	out_dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/bgc_hotspots/'

	cpg_limit = 0.001
	gc_limit = 0.01
	sp = 'LTF'

	# first match hotspots and coldspots
	# get the hotspots
	hotspot_file = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/spot2kb_flank40kb.seqldhot_validate_hotspots.csv'
	hot = get_hotspots(hotspot_file, sp)
	# get the coldspots
	coldspot_file = '/mnt/gluster/home/sonal.singhal1/%s/analysis/hotspots/%s.ldhelmet_unvalidated_coldspots.csv' % (sp, sp)
	cold = get_coldspots(coldspot_file)
	# then get cpg / CG data
	hot = get_cg_data(genome, hot)
	cold = get_cg_data(genome, cold)
		
	hot, cold = match_hotspots(hot, cold, cpg_limit, gc_limit)
	print_hotspots(out_dir, hot, cold)
	pickle.dump(hot, open('%shotspots.all_finches.pickle' % (out_dir), 'w'))
        pickle.dump(cold, open('%scoldspots.all_finches.pickle' % (out_dir), 'w'))


if __name__ == "__main__":
    main()
