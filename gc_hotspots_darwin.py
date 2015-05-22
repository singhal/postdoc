import re
import subprocess
from itertools import izip
import pandas as pd

def get_hs_cs(d):
	matches = {}
	
	for chr, hot_mid in zip(d.chr, d.hot_mid):
		if chr not in matches:
			matches[chr] = {}
		matches[chr][hot_mid] = 'hot'

	for chr, cold_mid in zip(d.chr, d.cold_mid):
		matches[chr][cold_mid] = 'cold'
	
	return matches


def get_ancestral_seq(genomes, matches, chr, spot):
	seq = {}

	for genome in genomes:
		start = spot - 5000
		end = spot + 5000
		seq[genome] = '' 

		seqcall = subprocess.Popen('samtools faidx %s %s:%s-%s' % (genomes[genome], chr, start, end), shell=True, stdout=subprocess.PIPE)
		for l in seqcall.stdout:
			if not re.match('>', l):
				seq[genome] += l.upper().rstrip()

	for genome, sequence in seq.items():
		seq[genome] = list(sequence)

	return seq


def get_lineage_mutations(seq, counts):
	# mutation types
	types = {       'A': {'A': 'AT_AT', 'G': 'AT_GC', 'T': 'AT_AT', 'C': 'AT_GC'},
                        'C': {'A': 'GC_AT', 'G': 'GC_GC', 'T': 'GC_AT', 'C': 'GC_GC'},
                        'T': {'A': 'AT_AT', 'G': 'AT_GC', 'T': 'AT_AT', 'C': 'AT_GC'},
                        'G': {'A': 'GC_AT', 'G': 'GC_GC', 'T': 'GC_AT', 'C': 'GC_GC'},
                }

	for ix, (anc, zf, ltf, dbf, gmag, gfor) in enumerate(izip(seq['anc'], seq['ZF'], seq['LTF'], seq['DBF'], seq['Gmag'], seq['Gfor'])):
		tuple1anc = tuple2anc = tuple1der = tuple2der = 'NN'
                if (ix - 1) >= 0:
                        tuple1anc = seq['anc'][ix - 1] + anc
			tuple1der = seq['Gfor'][ix - 1] + gfor
                if (ix + 1) < len(seq['anc']):
                        tuple2anc = anc + seq['anc'][ix + 1]
			tuple2der = gfor + seq['Gfor'][ix + 1]

		cpg = False
		for tuple in [tuple1anc, tuple1der, tuple2anc, tuple2der]:
			if tuple in ['CG', 'NG', 'CN']:
				cpg = True
		
		if anc != 'N' and gfor in ['A', 'T', 'C', 'G']:
			if not cpg:
				if anc in ['A', 'T']:
					counts[int(ix/100.)]['ancAT'] += 1
				elif anc in ['G', 'C']:
					counts[int(ix/100.)]['ancGC'] += 1
				if gfor != anc:
					if gfor not in [zf, ltf, dbf, gmag]:
						type = types[anc][gfor]
						if type in ['AT_GC', 'GC_AT']:
							counts[int(ix/100.)][type] += 1
	return counts


def main():
	genomes = {	'LTF': '/mnt/gluster/home/sonal.singhal1/reference/LTF_reference.fa',
			'ZF': '/mnt/gluster/home/sonal.singhal1/reference/ZF_reference.fa',
			'DBF': '/mnt/gluster/home/sonal.singhal1/reference/DBF_reference.fa',
			'anc': '/mnt/gluster/home/sonal.singhal1/reference/ancestral_genome.all_finches.fa',
			'Gfor': '/mnt/gluster/home/sonal.singhal1/Darwin/g_fortis/Geospiza_fortis.fasta',
			'Gmag': '/mnt/gluster/home/sonal.singhal1/Darwin/g_magnirostris/Geospiza_magnirostris.fasta'
		}
	
	out = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/bgc_hotspots/gfortis.hotspots_noCpG.csv'
	o = open(out, 'w')
	o.write('chr,midpoint,spot_type,bin,ancAT,ancGC,AT_GC,GC_AT\n')
	spots = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/bgc_hotspots/matched_hotspots.01GC.001CpG.all_finches.csv'
	spots = pd.read_csv(spots)
	
	# get matched hotspots and coldspots
	matches = get_hs_cs(spots)

	for chr in matches:
		for spot in matches[chr]:
			counts = {}
			for count in range(0,101):
				counts[count] = {'ancGC': 0, 'ancAT': 0, 'AT_GC': 0, 'GC_AT': 0}

			seq = get_ancestral_seq(genomes, matches, chr, spot)	
			# identify mutations in focal lineage
			# and get ancestral counts
			counts = get_lineage_mutations(seq, counts)	
			for bin in counts:
				o.write('%s,%s,%s,%s,%s,%s,%s,%s\n' % (chr, spot, matches[chr][spot], bin, counts[bin]['ancAT'], counts[bin]['ancGC'], counts[bin]['AT_GC'], counts[bin]['GC_AT']))
	o.close()
	

if __name__ == "__main__":
    main()
