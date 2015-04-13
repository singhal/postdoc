import subprocess
import pickle
import re
import random
import numpy

def get_sequence(sp):
	spots = pickle.load(open('/mnt/gluster/home/sonal.singhal1/ZF/analysis/motifs/%shotspots.pickle' % sp, "r"))
	hot = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/motifs/%s_hot_1kb_sequences.subset.fa' % (sp)
	cold = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/motifs/%s_cold_1kb_sequences.subset.fa' % (sp)
	h = open(hot, 'w')
	c = open(cold, 'w')
	
	num_spots = 0
	for chr in spots:
		for mid in spots[chr]:
			if spots[chr][mid]['match']:
				num_spots += 1
	num_samples = 1000

	for chr in spots:
		for mid in spots[chr]:
			if spots[chr][mid]['match']:
				mid = int(mid)
				start = mid - 500
				end = mid + 500

				match_mid = spots[chr][mid]['match']
				match_start = match_mid - 500
				match_end = match_mid + 500

				if random.random() < (num_samples / float(num_spots)):
					seq_call = subprocess.Popen('samtools faidx %s %s:%s-%s' % (genome, chr, start, end), shell=True, stdout=subprocess.PIPE)
					seq = '' 
					for l in seq_call.stdout:
						if not re.match('>', l):
							seq += l.rstrip().upper()
					h.write('>%s_%s\n%s\n' % (chr, mid, seq))

					seq_call = subprocess.Popen('samtools faidx %s %s:%s-%s' % (genome, chr, match_start, match_end), shell=True, stdout=subprocess.PIPE)
                                        seq = '' 
                                        for l in seq_call.stdout:
                                                if not re.match('>', l):
                                                        seq += l.rstrip().upper()
                                        c.write('>%s_%s\n%s\n' % (chr, match_mid, seq))
	h.close()
	c.close()

genome = '/mnt/gluster/home/sonal.singhal1/reference/ancestral_genome.fa'
get_sequence('ZF')
get_sequence('LTF')
