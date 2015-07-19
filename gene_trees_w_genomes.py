import re
import subprocess
import random
import numpy as np
import glob

def get_alignments(n_align, chr_lengths, genome):
	genome_length = sum(chr_lengths.values())
	chr_prop = {}
	loci = {}
	for chr, length in chr_lengths.items():
		chr_prop[ chr ] = round( n_align * length / float(genome_length))
	for chr, n_loci in chr_prop.items():
		tot_n = 0
		while tot_n < n_loci:
			start = random.randint(0, chr_lengths[chr])
			end = start + 999
			seq = subprocess.Popen('samtools faidx %s %s:%s-%s' % (genome, chr, start, end), shell=True, stdout=subprocess.PIPE)			
			n_count = 0
			for l in seq.stdout:
				if not re.match('>', l):
					n_count += list(l.rstrip()).count('0')

			if n_count > 750:
				tot_n += 1
				if chr not in loci:
					loci[chr] = dict()
				loci[chr][start] = 1	
	return loci

def make_fasta_files(loci, names, dir):
	blat_file = '%ssequences_to_blast.fa' % dir
	b_out = open(blat_file, 'w')

	for chr in loci:
		for start in loci[chr]:
			end = start + 999
			out_file = '%s%s_%s_%s.fasta' % (dir, chr, start, end)
			out_f = open(out_file, 'w')
			for species in names:
				genome = '/mnt/gluster/home/sonal.singhal1/gene_trees/chromosomes/%s_%s_haplotypes.fasta' % (species, chr)
				for haplo in names[species]:
					seq = subprocess.Popen('samtools faidx %s %s:%s-%s' % (genome, haplo, start, end), shell=True, stdout=subprocess.PIPE) 
					for l in seq.stdout:
                                		if re.match('>', l):
                                        		out_f.write('>%s\n' % names[species][haplo])
							if haplo == 'haplo0' and species == 'ZF':
								b_out.write('>%s_%s_%s\n' % (chr, start, end))
						else:
							out_f.write(l)
							if haplo == 'haplo0' and species == 'ZF':
								b_out.write(l)

	out_f.close()
	b_out.close()
	return blat_file


def blast_file_genome(dir, genomes, blast_file):
	results = {}
	for genome in genomes:
		out = '%s%s.blastn.out' % (dir, genome)
		results[genome] = out
		
		subprocess.call('blastn -query %s -db %s -out %s -max_target_seqs 1 -num_alignments 1 -outfmt 6 -num_threads 8' % (blast_file, genomes[genome], out), shell=True)
	return results


def parse_results(out_file):
	f = open(out_file, 'r')
	coords = {}
	for l in f:
		d = re.split('\s+', l.rstrip())
		if float(d[10]) < 1e-20:
			id = d[0]
			start = min(int(d[8]), int(d[9]))
			end = max(int(d[8]), int(d[9]))
			orient = '+'
			if int(d[8]) > int(d[9]):
				orient = '-'
			
			if id not in coords:
				coords[id] = {'contig': d[1], 'start': start, 'end': end, 'orient': orient}
			else:
				if d[1] == coords[id]['contig']:
					if start < coords[id]['start']:
						coords[id]['start'] = start
					if end > coords[id]['end']:
                                                coords[id]['end'] = end
	return coords


def rev_comp(seq):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
        return reverse_complement


def add_to_file(ids, coord, dir, genome, name):
	for id in coord:
		if id not in ids:
			ids[id] = 0
		ids[id] += 1
		
		fasta = '%s%s.fasta' % (dir, id)
		f = open(fasta, 'a')
		seq = subprocess.Popen('samtools faidx %s %s:%s-%s' % (genome, coord[id]['contig'], coord[id]['start'], coord[id]['end']), shell=True, stdout=subprocess.PIPE)
		s = ''
		for l in seq.stdout:
			if not re.match('>', l):
				s += l.rstrip()
		if coord[id]['orient'] == '-':
			s = rev_comp(s)
		f.write('>%s\n%s\n' % (name, s))
		f.close()
	return ids
		

def main():
	chr_lengths = {	'chr10': 20806668, 'chr11': 21403021, 'chr12': 21576510, 'chr13': 16962381,
			'chr14': 16419078, 'chr15': 14428146, 'chr1A': 73657157,
			'chr1': 118548696, 'chr2': 156412533, 'chr3': 112617285,
			'chr4A': 20704505, 'chr4': 69780378, 'chr5': 62374962, 'chr6': 36305782,
			'chr7': 39844632, 'chr8': 27993427, 'chr9': 27241186}

	names = {'LTF': {'haplo0': 'LTFh_73783a', 'haplo15': 'LTFh_73942b', 'haplo25': 'LTFa_G163b', 'haplo39': 'LTFa_W2994b'},
                'ZF': {'haplo0': 'ZF_26462a', 'haplo10': 'ZF_26792b', 'haplo25': 'ZF_28313b', 'haplo30': 'ZF_28402a'},
                'DBF': {'haplo0': 'DBFa', 'haplo1': 'DBFb'}}	

	genomes = {	'ficedula': '/mnt/gluster/home/sonal.singhal1/ficedula/Ficedula_albicollis.FicAlb15.fa',
			'gfortis': '/mnt/gluster/home/sonal.singhal1/Darwin/g_fortis/geoFor1.masked.fa' }		

	dir = '/mnt/gluster/home/sonal.singhal1/gene_trees2/'
	genome = '/mnt/gluster/home/sonal.singhal1/reference/all.masked_genome.fa'
	loci = get_alignments(1200, chr_lengths, genome) 
	blast_file = make_fasta_files(loci, names, dir)
	results = blast_file_genome(dir, genomes, blast_file)
	fic_coord = parse_results(results['ficedula'])
	geo_coord = parse_results(results['gfortis'])
	ids = {}
	ids = add_to_file(ids, fic_coord, dir, genomes['ficedula'], 'ficedula')
	ids = add_to_file(ids, geo_coord, dir, genomes['gfortis'], 'gfortis')
	for id in ids:
		if ids[id] == 2:
			print 'muscle -in %s%s.fasta -out %s%s.fasta.aln' % (dir, id, dir, id)

if __name__ == "__main__":
    main()


